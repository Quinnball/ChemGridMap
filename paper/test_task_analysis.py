import json

import pytest

from analyze_tasks import analyze, score


def test_empty_responses_do_not_become_successful_user_study(tmp_path):
    key=tmp_path/'keys.json'; key.write_text('{}')
    response=tmp_path/'responses.json'; response.write_text(json.dumps({'schema':1,'observations':[]}))
    with pytest.raises(ValueError,match='No human observations'):
        analyze([response],key,tmp_path/'output')
    assert not (tmp_path/'output').exists()


def test_omitted_molecular_answers_are_scored_as_missing_not_removed():
    key={'A':{'contains_aromatic_nitrogen':True,'activity_class':'active'},
         'B':{'contains_aromatic_nitrogen':False,'activity_class':'inactive'}}
    observation={'kind':'structure','answers':[{'id':'A','motif':'Yes','activity_class':'active'}]}
    assert score(observation,key) == .5


def test_blank_record_answers_do_not_count_as_zero():
    key={'n_records':2,'median_pchembl':7.185,'n_documents':1,'is_conflicted':True}
    observation={'kind':'records','answers':{'n_records':'2','median_pchembl':'7.185','n_documents':'','is_conflicted':''}}
    assert score(observation,key) == .5
