import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from finemap_qc import annotate, cli, io, stats, uncertainty
from finemap_qc.references import high_pip_variants
from test_finemap_qc import GTF


def read(tmp_path, text, **kw):
    p = tmp_path / 'cs.tsv'
    p.write_text(text)
    return io.read_credible_sets(p, **kw)


@pytest.mark.parametrize('pip', ['-0.1', '1.1', 'nan', 'inf', '-inf', 'bad', ''])
def test_invalid_pip(tmp_path, pip):
    with pytest.raises(ValueError):
        read(tmp_path, f'cs_id\tvariant\tpip\nA\t1_10_A_G\t{pip}\n')
    with pytest.raises(ValueError):
        read(tmp_path, f'A\t1_10_A_G\t{pip}\n')


@pytest.mark.parametrize('text', ['', 'a\tb\n', 'cs_id\tvariant\tpip\n',
                                  'cs_id\tvariant\tpip\n\t1_10_A_G\t0.9\n'])
def test_empty_or_malformed(tmp_path, text):
    with pytest.raises(ValueError):
        read(tmp_path, text)


def test_duplicates_and_traits(tmp_path):
    text = 'cs_id\tvariant\tpip\tmolecular_trait_id\nA\t1_10_A_G\t0.9\tG1\nA\tchr1:10:A:G\t0.9\tG1\n'
    with pytest.raises(ValueError, match='duplicate'):
        read(tmp_path, text)
    cs = read(tmp_path, text, duplicate_policy='collapse-identical')
    assert len(cs) == 1 and cs.attrs['duplicate_rows_removed'] == 1
    with pytest.raises(ValueError, match='duplicate'):
        read(tmp_path, text.replace('chr1:10:A:G\t0.9', 'chr1:10:A:G\t0.8'), duplicate_policy='collapse-identical')
    with pytest.raises(ValueError, match='inconsistent trait'):
        read(tmp_path, text.replace('chr1:10:A:G', '1_20_A_G').replace('0.9\tG1\n', '0.9\tG2\n', 1))


@pytest.mark.parametrize('variant', ['1_0_A_G', '1_1.5_A_G', '1_inf_A_G', '1_1e3_A_G', '1_10__G'])
def test_invalid_position_or_allele(tmp_path, variant):
    with pytest.raises(ValueError, match='unparsable'):
        read(tmp_path, f'A\t{variant}\t0.5\n')


def test_named_catalogue_order_and_boundaries(tmp_path):
    cs = read(tmp_path, 'molecular_trait_id\tcs_id\tvariant\tpip\tlocus_id\nG\tA\t1_10_A_G\t0\tL\nG\tA\t1_20_A_G\t1\tL\n')
    assert cs.trait.tolist() == ['G', 'G']
    assert cs.locus_id.tolist() == ['L', 'L']
    assert io.normalize_variant_ids(['1_001_a_g']).variant.iloc[0] == 'chr1_1_A_G'


def test_canonical_splice_and_single_exon(tmp_path):
    p = tmp_path / 'a.gtf'
    p.write_text(GTF)
    g = annotate.load_gene_features(p, True)
    assert g['splice']['1'].tolist() == [1200, 3000, 3100, 4800]
    # Select the minus-strand transcript, then make it single-exon.
    p2 = tmp_path / 'minus.gtf'
    minus = GTF.replace('tag "Ensembl_canonical";', '').replace('transcript_id "T2";', 'transcript_id "T2"; tag "Ensembl_canonical";')
    p2.write_text(minus)
    assert annotate.load_gene_features(p2, True)['splice']['1'].tolist() == [21000, 29000]
    p3 = tmp_path / 'single.gtf'
    p3.write_text('\n'.join(line for line in minus.splitlines() if '\texon\t29000' not in line))
    assert annotate.load_gene_features(p3, True)['splice'] == {}


def test_cluster_components():
    cs = pd.DataFrame({'cs_id':['A','A','B','B','C','D'], 'variant':['v1','v2','v2','v3','v4','v5']})
    c = uncertainty.variant_clusters(cs)
    assert c['v1'] == c['v3'] and c['v1'] != c['v4']
    cs['locus_id'] = ['L1','L1','L2','L2','L2','L3']
    c = uncertainty.variant_clusters(cs)
    assert c['v1'] == c['v4'] and c['v4'] != c['v5']
    pd.testing.assert_series_equal(c.sort_index(), uncertainty.variant_clusters(cs.sample(frac=1, random_state=9)).sort_index())


def test_bootstrap_known_answer_and_missing():
    ids = ['h1','l1','h2','l2','h3','l3']
    annot = pd.DataFrame({'b':[1,1,1,1,1,1], 'q':[3,1,3,1,3,1]}, index=ids)
    groups = {'high':set(ids[::2]), 'low_all':set(ids[1::2]), 'low_anchored':set(ids[1::2])}
    clusters = pd.Series(['A','A','B','B','C','C'],index=ids)
    b = uncertainty.cluster_intervals(stats.binary_enrichment(annot,groups,['b']),annot,groups,clusters,True,100)
    assert b.fe_ci_low.tolist() == [1,1] and b.fe_ci_high.tolist() == [1,1]
    q = uncertainty.cluster_intervals(stats.quantitative_comparison(annot,groups,['q']),annot,groups,clusters,False,100)
    assert q.median_diff_ci_low.tolist() == [2,2] and q.median_diff_ci_high.tolist() == [2,2]
    pd.testing.assert_frame_equal(q, uncertainty.cluster_intervals(stats.quantitative_comparison(annot,groups,['q']),annot,groups,clusters,False,100))
    annot.loc['l1','b'] = np.nan
    assert stats.binary_enrichment(annot,groups,['b']).n_low.tolist() == [2,2]
    clusters[:] = 'A'
    b = uncertainty.cluster_intervals(b,annot,groups,clusters,True,100)
    assert b.fe_ci_low.isna().all() and (b.ci_status == 'insufficient_clusters').all()


def test_bootstrap_undefined_denominator():
    annot = pd.DataFrame({'b':[1,0,1,0]},index=['h1','l1','h2','l2'])
    groups={'high':{'h1','h2'},'low_all':{'l1','l2'},'low_anchored':{'l1','l2'}}
    c = pd.Series(['A','A','B','B'],index=annot.index)
    b = uncertainty.cluster_intervals(stats.binary_enrichment(annot,groups,['b']),annot,groups,c,True,100)
    assert (b.ci_status == 'undefined_replicates').all() and b.fe_ci_low.isna().all()


def test_weighted_median_matches_expansion():
    values = np.array([1,2,10,25])
    for counts in ([0,2,1,0], [1,0,0,1], [3,2,0,1]):
        assert uncertainty._weighted_median(values, counts) == np.median(np.repeat(values,counts))


@pytest.mark.parametrize('high,low,anchor', [(0.8,0.9,0.5),(np.nan,0.05,0.5),(0.8,0.05,2)])
def test_bad_thresholds(high,low,anchor):
    with pytest.raises(ValueError):
        stats.define_groups(pd.DataFrame(),high,low,anchor)


def load_script(name):
    spec = importlib.util.spec_from_file_location(name,Path(__file__).parents[1]/'scripts'/f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_reference_study_pooling_and_exclusion(tmp_path):
    paths=[]
    for i in range(3):
        p=tmp_path/f'{i}.tsv'
        p.write_text('variant\tpip\n1_10_A_G\t0.9\n1_20_A_G\t0.8\n')
        paths.append(str(p))
    out=tmp_path/'rep.txt'
    builder=load_script('build_replicating_set')
    builder.main(paths[:2]+['--study-ids','A','A','--out',str(out)])
    assert out.read_text()==''
    builder.main(paths+['--study-ids','A','A','B','--out',str(out)])
    assert out.read_text()=='chr1_10_A_G\n'
    builder.main(paths+['--study-ids','A','A','B','--exclude-study','B','--out',str(out)])
    assert out.read_text()==''


def test_reference_invalid_ids_and_pips():
    for v,p in [('rs123',0.9),('1_10_A_G',1.1),('1_10_A_G',np.nan)]:
        with pytest.raises(ValueError):
            high_pip_variants(pd.DataFrame({'variant':[v],'pip':[p]}),0.8)


def test_real_catalogue_fixture():
    path = Path(__file__).parent/'data'/'QTD000001.one_cs.tsv'
    raw = pd.read_csv(path,sep='\t')
    with pytest.raises(ValueError,match='duplicate'):
        io.read_credible_sets(path)
    cs = io.read_credible_sets(path,'collapse-identical')
    assert len(cs) == int(raw.cs_size.iloc[0]) == 53
    assert cs.cs_id.nunique() == 1 and cs.trait.nunique() == 1
    assert cs.attrs['duplicate_rows_removed'] == len(raw)-53


def test_catalogue_reference_cache_and_provenance(tmp_path):
    import json
    builder = load_script('download_references')
    rows=[]
    for i,study in enumerate(['A','A','B']):
        source=tmp_path/f'source{i}.tsv'
        source.write_text('variant\tpip\n1_10_A_G\t0.9\n1_20_A_G\t0.8\n')
        rows.append(dict(dataset_id=f'D{i}', study_id=study,quant_method='ge',
                         study_label=study,sample_group='tissue',ftp_cs_path=str(source)))
    meta=tmp_path/'meta.tsv'
    pd.DataFrame(rows).to_csv(meta,sep='\t',index=False)
    out=builder.eqtl_catalogue(tmp_path,str(meta),['ge'],0.8,2,set())
    assert (tmp_path/out['eqtlcat_ge_rep']).read_text()=='chr1_10_A_G\n'
    provenance=json.loads((tmp_path/'eqtlcat_ge_rep.provenance.json').read_text())
    assert len(provenance['datasets'])==3 and provenance['min_studies']==2
    # Cached rows suffice after source files are removed.
    for row in rows:
        Path(row['ftp_cs_path']).unlink()
    builder.eqtl_catalogue(tmp_path,str(meta),['ge'],0.8,2,{'B'})
    assert (tmp_path/'eqtlcat_ge_rep.txt').read_text()==''


def test_cli_cluster_outputs_and_group_priority(tmp_path):
    p=tmp_path/'cs.tsv'
    rows=[]
    for k in range(3):
        rows += [(f'A{k}',f'1_{100+k*10}_A_G',.9),(f'A{k}',f'1_{101+k*10}_A_G',.01)]
    pd.DataFrame(rows,columns=['cs_id','variant','pip']).to_csv(p,sep='\t',index=False)
    bed=tmp_path/'a.bed'; bed.write_text('chr1\t0\t1000\n')
    out=tmp_path/'out'
    cli.main(['--credible-sets',str(p),'--out',str(out),'--bed',f'open={bed}','--bootstrap-replicates','100'])
    tab=pd.read_csv(out/'binary_enrichment.tsv',sep='\t')
    assert tab.ci_method.tolist()==['cluster_percentile']*2
    assert tab.fe_ci_low.tolist()==[1,1]
    ann=pd.read_csv(out/'variant_annotations.tsv.gz',sep='\t')
    assert (ann.group=='low_anchored').sum()==3
    assert 'cluster bootstrap' in (out/'report.html').read_text()


def test_cluster_bootstrap_matches_explicit_resampling():
    ids=['h1','l1','l2','h2','h3','l3']
    a=pd.DataFrame({'q':[1,5,9,20,30,11]},index=ids)
    g={'high':{'h1','h2','h3'},'low_all':{'l1','l2','l3'},'low_anchored':{'l1','l2','l3'}}
    clusters=pd.Series(['A','A','A','B','C','C'],index=ids)
    result=uncertainty.cluster_intervals(stats.quantitative_comparison(a,g,['q']),a,g,clusters,False,100,42)
    rng=np.random.default_rng(42); expected=[]
    for _ in range(100):
        draw=rng.integers(3,size=3)
        h=[];l=[]
        for k in draw:
            for v in clusters.index[clusters==['A','B','C'][k]]:
                (h if v in g['high'] else l).append(a.loc[v,'q'])
        expected.append(np.median(h)-np.median(l) if h and l else np.nan)
    assert result.bootstrap_valid.iloc[0]==np.isfinite(expected).sum()
    if np.isfinite(expected).all():
        np.testing.assert_allclose(result.loc[0,['median_diff_ci_low','median_diff_ci_high']].astype(float),np.quantile(expected,[.025,.975]))
    else:
        assert result.ci_status.iloc[0]=='undefined_replicates'


def test_strict_group_threshold_boundaries():
    cs=pd.DataFrame({'cs_id':['A','A','A','B','B','C','C'],
                     'variant':['h','at_high','at_low','at_anchor','b_low','above_anchor','c_low'],
                     'pip':[.9,.8,.05,.5,.01,.50001,.01]})
    groups=stats.define_groups(cs)
    assert groups['high']=={'h'}
    assert groups['low_all']=={'b_low','c_low'}
    assert groups['low_anchored']=={'c_low'}


def test_inconsistent_or_missing_locus(tmp_path):
    for locus in ['L2','']:
        with pytest.raises(ValueError,match='locus_id'):
            read(tmp_path, f'cs_id\tvariant\tpip\tlocus_id\nA\t1_10_A_G\t0.9\tL1\nA\t1_20_A_G\t0.01\t{locus}\n')
