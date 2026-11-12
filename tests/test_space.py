"""Optional local Space integration test: install gradio==6.31.0 to run."""
import importlib.util
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest


def test_space_duplicate_option_and_report(tmp_path, monkeypatch):
    pytest.importorskip('gradio')
    from fastapi.testclient import TestClient
    root = Path(__file__).parents[1]
    stage = tmp_path/'stage'; stage.mkdir()
    shutil.copy(root/'space/app.py', stage/'app.py')
    shutil.copy(root/'space/index.html', stage/'index.html')
    shutil.copytree(root/'examples', stage/'examples')
    monkeypatch.setenv('FINEMAP_QC_REFS_REPO','')
    monkeypatch.setenv('GRADIO_ANALYTICS_ENABLED','False')
    spec = importlib.util.spec_from_file_location('space_test_app', stage/'app.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    module.WORK = tmp_path/'runs'; module.WORK.mkdir()
    # Keep the real analysis, but execute it synchronously for deterministic API assertions.
    monkeypatch.setattr(module.threading, 'Thread', lambda target,args,daemon: SimpleNamespace(start=lambda:target(*args)))
    client = TestClient(module.app)
    assert client.get('/').status_code == 200
    fixture=(root/'tests/data/QTD000001.one_cs.tsv').read_bytes()
    for collapse, expected in [('false','error'),('true','done'),(None,'done')]:  # None: default merges
        data={'use_refs':'false'} | ({} if collapse is None else {'collapse_duplicates':collapse})
        response=client.post('/api/run', data=data,
                             files={'credible_sets':('cs.tsv',fixture,'text/tab-separated-values')})
        assert response.status_code==200
        rid=response.json()['id']
        assert client.get(f'/api/runs/{rid}').json()['status']==expected
        if expected=='done':
            page=client.get(f'/runs/{rid}/report.html')
            assert page.status_code==200 and 'collapse-identical' in page.text
            assert client.get(f'/runs/{rid}/results.zip').status_code==200
