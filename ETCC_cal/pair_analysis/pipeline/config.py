"""パイプラインの設定JSONを読み、カットの組（cut set）を展開し、各ステージの保存先を決める。

設定JSONの中身は config/ の例と README.md を参照。相対パスはリポジトリ直下（analysis_sm3_FM）が基準。
"""
import copy
import json
from pathlib import Path

from settings import load_config as load_etcc_config
from settings import repo_path

STAGE_DIRS = {
    0: "stage0_inputs",
    1: "stage1_vertices",
    2: "stage2_branch_fit",
    3: "stage3_classify",
    4: "stage4_modulation",
    5: "stage5_event_display",
}


def deep_merge(base, override):
    """辞書を入れ子ごと上書きした新しい辞書を返す（base は変えない）。"""
    merged = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def strip_comments(value):
    """"_" で始まるキー（説明用）を取り除く。"""
    if isinstance(value, dict):
        return {k: strip_comments(v) for k, v in value.items() if not k.startswith("_")}
    if isinstance(value, list):
        return [strip_comments(v) for v in value]
    return value


class Pipeline:
    """設定JSON 1つ分の情報。各ステージはこれだけを受け取る。"""

    def __init__(self, path):
        self.path = Path(path).resolve()
        self.raw = json.loads(self.path.read_text())
        settings = strip_comments(self.raw)
        self.name = settings["name"]
        self.output_dir = repo_path(settings["output_dir"])
        self.beam_axis_file = repo_path(settings["beam_axis_fit"])
        self.datasets = settings["datasets"]
        self.scan = settings["scan"]
        self.analysis = settings["analysis"]
        self.comparisons = settings["comparisons"]
        self.event_display = settings.get("event_display", {"enabled": False})
        self.samples = self._build_samples(settings.get("samples", {}))
        self.cut_sets = self._expand_cut_sets(settings["cut_sets"])
        # Z 側除外を外す cut set が1つでもあれば、stage2 で除外にかかったイベントも読む（読む量が約3倍）
        self.scan_includes_z_vetoed = any(not cuts["z_veto"] for cuts in self.cut_sets.values())
        self._check()

    # ---- 設定の展開 ----
    def _build_samples(self, extra):
        """データセット（全ラン）と、追加で定義したランの部分集合を、同じ形 {dataset, runs, label} にそろえる。"""
        samples = {name: {"dataset": name, "runs": list(d["runs"]), "label": d.get("label", name)}
                   for name, d in self.datasets.items()}
        for name, sample in extra.items():
            dataset = self.datasets[sample["dataset"]]
            runs = sample.get("runs", dataset["runs"])
            missing = set(runs) - set(dataset["runs"])
            if missing:
                raise ValueError(f"sample {name}: runs {sorted(missing)} are not in dataset {sample['dataset']}")
            samples[name] = {"dataset": sample["dataset"], "runs": list(runs), "label": sample.get("label", name)}
        return samples

    @staticmethod
    def _expand_cut_sets(definitions):
        """"base" で別の cut set を継承し、差分だけ書けるようにする。"""
        expanded = {}

        def resolve(name, chain=()):
            if name in expanded:
                return expanded[name]
            if name in chain:
                raise ValueError(f"cut set inheritance loop: {' -> '.join(chain + (name,))}")
            definition = dict(definitions[name])
            base = definition.pop("base", None)
            cuts = deep_merge(resolve(base, chain + (name,)), definition) if base else definition
            expanded[name] = cuts
            return cuts

        for name in definitions:
            resolve(name)
        return expanded

    def _check(self):
        """スキャンした範囲より広いカットを指定していないか確かめる。"""
        for name, cuts in self.cut_sets.items():
            if cuts["beam_radius_mm"] > self.scan["beam_radius_mm"]:
                raise ValueError(f"cut set {name}: beam_radius_mm {cuts['beam_radius_mm']} > scan.beam_radius_mm {self.scan['beam_radius_mm']}")
        for comparison in self.comparisons:
            for key in ("numerator", "denominator"):
                if comparison[key] not in self.samples:
                    raise ValueError(f"comparison {comparison['name']}: unknown sample {comparison[key]}")

    # ---- 入出力の場所 ----
    def stage_dir(self, stage, *parts):
        """<output_dir>/stageN_xxx/parts... を作って返す。"""
        path = self.output_dir / STAGE_DIRS[stage]
        for part in parts:
            path = path / part
        path.mkdir(parents=True, exist_ok=True)
        return path

    def etcc_config(self, dataset):
        """raw の場所・較正・散乱点CSV・Z側除外の場所を持つ、既存の ETCC_cal 設定。"""
        return load_etcc_config(repo_path(self.datasets[dataset]["etcc_config"]))

    def selected_cut_sets(self, names=None):
        names = names or list(self.cut_sets)
        unknown = set(names) - set(self.cut_sets)
        if unknown:
            raise ValueError(f"unknown cut sets: {sorted(unknown)}")
        return {name: self.cut_sets[name] for name in names}
