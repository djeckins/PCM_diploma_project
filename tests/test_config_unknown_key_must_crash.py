"""Defect class: a typo in a config key silently falls back to the default."""

import pytest
import yaml

from pcm2 import config


def _minimal(tmp_path, extra=None, origins=None):
    doc = {
        "system": {"id": "t"},
        "data": {"data_root": "/x", "conditions": [
            {"id": "c1", "replicas": [
                {"id": "r1", "topology": "/t.tpr", "trajectory": "/t.xtc",
                 "lineage": "c1/r1"}]}]},
        "labels": {"tau_ps": [1000.0], "primary_tau_ps": 1000.0},
        "model": {"seed": 1},
    }
    if extra:
        doc.update(extra)
    if origins:
        doc["origins"] = origins
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(doc, allow_unicode=True))
    return p


def test_unknown_key_crashes(tmp_path):
    p = _minimal(tmp_path, extra={"labelz": {"tau_ps": [1.0]}})
    with pytest.raises(config.ConfigError, match="unknown keys"):
        config.load(p)


def test_unknown_nested_key_crashes(tmp_path):
    doc = yaml.safe_load(_minimal(tmp_path).read_text())
    doc["labels"]["tau_pss"] = [1.0]
    p = tmp_path / "cfg2.yaml"
    p.write_text(yaml.safe_dump(doc, allow_unicode=True))
    with pytest.raises(config.ConfigError, match="unknown keys"):
        config.load(p)


def test_valid_config_loads(tmp_path):
    cfg = config.load(_minimal(tmp_path))
    assert cfg["labels.primary_tau_ps"] == 1000.0
    # The three provenances are distinguishable: an untouched key is default.
    assert cfg.origin_of("data.stride")[0] == "default"


def test_origin_without_basis_crashes(tmp_path):
    p = _minimal(tmp_path, origins={"model.seed": {"origin": "declared", "basis": ""}})
    with pytest.raises(config.ConfigError, match="empty basis"):
        config.load(p)


def test_primary_tau_must_be_in_sweep(tmp_path):
    doc = yaml.safe_load(_minimal(tmp_path).read_text())
    doc["labels"]["primary_tau_ps"] = 999.0
    p = tmp_path / "cfg3.yaml"
    p.write_text(yaml.safe_dump(doc, allow_unicode=True))
    with pytest.raises(config.ConfigError, match="primary_tau_ps"):
        config.load(p)


def test_detected_write_records_origin(tmp_path):
    cfg = config.load(_minimal(tmp_path))
    cfg.set_key("data.pore.axis", "z", "detected", "test: measured")
    assert cfg.origin_of("data.pore.axis") == ("detected", "test: measured")
    config.save(cfg)
    cfg2 = config.load(cfg.source_path)
    assert cfg2["data.pore.axis"] == "z"
    assert cfg2.origin_of("data.pore.axis")[0] == "detected"


def test_saved_file_does_not_freeze_defaults(tmp_path):
    """Defect class: a default value materialized into the file freezes and
    stops tracking the schema -- a later schema edit silently fails to apply."""
    import yaml as _yaml
    cfg = config.load(_minimal(tmp_path))
    config.save(cfg)
    doc = _yaml.safe_load(cfg.source_path.read_text())
    # A key holding the untouched schema value is absent from the file...
    assert "trees" not in doc.get("model", {}), \
        "default keys must not be materialized into the file"
    # ...yet on load it exists and equals the current schema value.
    cfg2 = config.load(cfg.source_path)
    assert cfg2["model.trees.depth_grid"] == \
        config.SCHEMA["model"]["trees"]["depth_grid"].default
    # Required keys, in contrast, are written.
    assert doc["system"]["id"] == "t"
