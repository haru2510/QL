"""図の文字を日本語にするための設定（README の説明用の図で使う）。"""
import matplotlib

JAPANESE_FONTS = ["Noto Sans CJK JP", "Droid Sans Japanese", "DejaVu Sans"]


def use_japanese_font():
    """日本語のフォントを使う。マイナス記号はフォントに無いことがあるので ASCII の '-' にする。"""
    matplotlib.rcParams["font.family"] = "sans-serif"
    matplotlib.rcParams["font.sans-serif"] = JAPANESE_FONTS + list(matplotlib.rcParams["font.sans-serif"])
    matplotlib.rcParams["axes.unicode_minus"] = False
