from __future__ import annotations

import math
import shutil
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

OUT_PATH = Path(
    "/Users/josephine/Research/Space-FALCON/_Proposals/[FY27]SURP/AERO-08_Preliminary_Results/plots/body_schematic.png"
)
VAULT_COPY = Path(
    "/Users/josephine/Research/Space-FALCON/_Knowledge/Ideas/Details/assets/aero-08-prelim-figures/body_schematic.png"
)

W, H = 2200, 1320
BG = (248, 251, 255)
CARD = (255, 255, 255)
OUTLINE = (214, 226, 242)
NAVY = (20, 33, 56)
TEXT = (44, 57, 78)
MUTED = (105, 117, 136)
BLUE = (49, 102, 255)
TEAL = (63, 172, 180)
LIGHT_BLUE = (232, 240, 252)
BODY = (199, 214, 236)
BODY_DARK = (161, 182, 211)
PANEL = (245, 157, 35)
PANEL_FILL = (255, 227, 178)
PANEL_SHADE = (241, 189, 97)
GRID = (214, 223, 238)


def font(size: int, bold: bool = False):
    candidates = [
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
        if bold
        else "/System/Library/Fonts/Supplemental/Arial.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
        "/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"
        if bold
        else "/System/Library/Fonts/Supplemental/Tahoma.ttf",
    ]
    for path in candidates:
        try:
            return ImageFont.truetype(path, size=size)
        except Exception:
            continue
    return ImageFont.load_default()


def text(draw: ImageDraw.ImageDraw, xy, s, *, size=34, fill=TEXT, bold=False, anchor="la", spacing=6):
    draw.multiline_text(
        xy,
        s,
        font=font(size, bold=bold),
        fill=fill,
        anchor=anchor,
        spacing=spacing,
    )


def arrow(draw, p1, p2, *, fill=BLUE, width=6, head=20):
    draw.line([p1, p2], fill=fill, width=width)
    ang = math.atan2(p2[1] - p1[1], p2[0] - p1[0])
    left = (p2[0] - head * math.cos(ang - 0.45), p2[1] - head * math.sin(ang - 0.45))
    right = (p2[0] - head * math.cos(ang + 0.45), p2[1] - head * math.sin(ang + 0.45))
    draw.polygon([p2, left, right], fill=fill)


def dimension(draw, p1, p2, label, *, fill=BLUE, size=24, offset=(0, 0), anchor="mm"):
    draw.line([p1, p2], fill=fill, width=5)
    tick = 18
    if p1[0] == p2[0]:
        draw.line([(p1[0] - tick, p1[1]), (p1[0] + tick, p1[1])], fill=fill, width=5)
        draw.line([(p2[0] - tick, p2[1]), (p2[0] + tick, p2[1])], fill=fill, width=5)
    else:
        draw.line([(p1[0], p1[1] - tick), (p1[0], p1[1] + tick)], fill=fill, width=5)
        draw.line([(p2[0], p2[1] - tick), (p2[0], p2[1] + tick)], fill=fill, width=5)
    cx = (p1[0] + p2[0]) / 2 + offset[0]
    cy = (p1[1] + p2[1]) / 2 + offset[1]
    text(draw, (cx, cy), label, size=size, fill=fill, anchor=anchor)


def rounded_box(draw, box, *, title: str | None = None):
    draw.rounded_rectangle(box, radius=30, fill=CARD, outline=OUTLINE, width=3)
    if title:
        x0, y0, _, _ = box
        text(draw, (x0 + 28, y0 + 24), title, size=30, bold=True, fill=NAVY)


def body_side_polygon(nose_x: int, cy: int, shoulder_x: int, aft_x: int, nose_r: int, aft_half: int):
    return [
        (nose_x, cy - int(0.46 * nose_r)),
        (shoulder_x, cy - int(0.86 * aft_half)),
        (aft_x, cy - aft_half),
        (aft_x + 24, cy - aft_half),
        (aft_x + 24, cy + aft_half),
        (aft_x, cy + aft_half),
        (shoulder_x, cy + int(0.86 * aft_half)),
        (nose_x, cy + int(0.46 * nose_r)),
    ]


def draw_side_body(draw: ImageDraw.ImageDraw, *, nose_x: int, cy: int, shoulder_x: int, aft_x: int, nose_r: int, aft_half: int):
    draw.polygon(
        body_side_polygon(nose_x, cy, shoulder_x, aft_x, nose_r, aft_half),
        fill=BODY,
        outline=NAVY,
    )
    draw.ellipse((nose_x - nose_r, cy - nose_r, nose_x + nose_r, cy + nose_r), fill=BODY, outline=NAVY, width=5)


def draw_state_card(draw: ImageDraw.ImageDraw, box, *, title: str, subtitle: str, deployed: bool):
    rounded_box(draw, box)
    x0, y0, x1, y1 = box
    text(draw, (x0 + 26, y0 + 24), title, size=32, bold=True, fill=NAVY)
    text(draw, (x0 + 26, y0 + 68), subtitle, size=22, fill=MUTED)

    cy = y0 + 265
    nose_x = x0 + 165
    shoulder_x = x0 + 335
    aft_x = x1 - 150
    nose_r = 62
    aft_half = 132

    draw.line([(x0 + 40, cy), (x1 - 40, cy)], fill=GRID, width=3)
    draw_side_body(draw, nose_x=nose_x, cy=cy, shoulder_x=shoulder_x, aft_x=aft_x, nose_r=nose_r, aft_half=aft_half)

    if deployed:
        upper_hinge = (aft_x - 235, cy - 112)
        lower_hinge = (aft_x - 235, cy + 112)
        upper_panel = [(upper_hinge[0], upper_hinge[1]), (aft_x + 40, cy - 205), (aft_x + 135, cy - 166), (upper_hinge[0] + 86, upper_hinge[1] + 16)]
        lower_panel = [(lower_hinge[0], lower_hinge[1]), (aft_x + 40, cy + 205), (aft_x + 135, cy + 166), (lower_hinge[0] + 86, lower_hinge[1] - 16)]
        draw.polygon(upper_panel, fill=PANEL_FILL, outline=PANEL)
        draw.polygon(lower_panel, fill=PANEL_FILL, outline=PANEL)
        draw.line([upper_hinge, (upper_hinge[0] + 34, upper_hinge[1] + 8)], fill=NAVY, width=8)
        draw.line([lower_hinge, (lower_hinge[0] + 34, lower_hinge[1] - 8)], fill=NAVY, width=8)
        arrow(draw, (x0 + 285, y0 + 138), (x0 + 438, y0 + 138), fill=TEAL, width=6, head=18)
        text(draw, (x0 + 285, y0 + 106), "Entry direction", size=21, bold=True, fill=TEAL)
        text(draw, (x0 + 28, y1 - 52), "2 symmetric broadside drag panels", size=23, bold=True, fill=PANEL)
    else:
        arrow(draw, (x0 + 285, y0 + 138), (x0 + 438, y0 + 138), fill=TEAL, width=6, head=18)
        text(draw, (x0 + 285, y0 + 106), "Entry direction", size=21, bold=True, fill=TEAL)
        text(draw, (x0 + 28, y1 - 52), "Passive body-only entry state", size=23, bold=True, fill=BLUE)


def draw_top_view(draw: ImageDraw.ImageDraw, box):
    rounded_box(draw, box, title="Top view: two symmetric panels")
    x0, y0, x1, y1 = box
    text(draw, (x0 + 28, y0 + 68), "Symmetric deployment adds drag, not lateral force.", size=21, fill=MUTED)

    cx = x0 + 330
    cy = y0 + 210
    nose_x = x0 + 120
    aft_x = x1 - 120
    half_front = 62
    half_back = 134

    draw.line([(x0 + 42, cy), (x1 - 42, cy)], fill=GRID, width=3)
    body = [
        (nose_x, cy - 40),
        (cx - 120, cy - half_front),
        (aft_x, cy - half_back),
        (aft_x + 22, cy - half_back),
        (aft_x + 22, cy + half_back),
        (aft_x, cy + half_back),
        (cx - 120, cy + half_front),
        (nose_x, cy + 40),
    ]
    draw.polygon(body, fill=BODY, outline=NAVY)
    draw.ellipse((nose_x - 54, cy - 54, nose_x + 54, cy + 54), fill=BODY, outline=NAVY, width=4)

    left_panel = [(cx + 18, cy - 116), (aft_x - 46, cy - 190), (aft_x + 16, cy - 172), (cx + 64, cy - 100)]
    right_panel = [(cx + 18, cy + 116), (aft_x - 46, cy + 190), (aft_x + 16, cy + 172), (cx + 64, cy + 100)]
    draw.polygon(left_panel, fill=PANEL_FILL, outline=PANEL)
    draw.polygon(right_panel, fill=PANEL_FILL, outline=PANEL)

    arrow(draw, (x0 + 64, y0 + 124), (x0 + 220, y0 + 124), fill=TEAL, width=6, head=18)
    text(draw, (x0 + 66, y0 + 92), "Entry direction", size=21, bold=True, fill=TEAL)
    text(draw, (x1 - 34, cy - 132), "Left panel", size=19, fill=PANEL, anchor="ra")
    text(draw, (x1 - 34, cy + 112), "Right panel", size=19, fill=PANEL, anchor="ra")


def draw_side_geometry(draw: ImageDraw.ImageDraw, box):
    rounded_box(draw, box, title="Analysis geometry used in the study")
    x0, y0, x1, y1 = box
    text(draw, (x0 + 28, y0 + 68), "SHIELD-informed 70° sphere-cone surrogate", size=21, fill=MUTED)

    cy = y0 + 220
    nose_x = x0 + 120
    shoulder_x = x0 + 280
    aft_x = x1 - 160
    nose_r = 52
    aft_half = 132

    draw.line([(x0 + 50, cy), (x1 - 50, cy)], fill=GRID, width=3)
    draw_side_body(draw, nose_x=nose_x, cy=cy, shoulder_x=shoulder_x, aft_x=aft_x, nose_r=nose_r, aft_half=aft_half)

    dimension(draw, (aft_x + 78, cy - aft_half), (aft_x + 78, cy + aft_half), "Base radius\n0.25 m", fill=BLUE, size=21, offset=(46, 0), anchor="lm")
    dimension(draw, (nose_x - nose_r, cy + 112), (nose_x + nose_r, cy + 112), "Nose radius 0.08 m", fill=BLUE, size=21, offset=(0, 36))
    arrow(draw, (x1 - 220, y0 + 108), (x0 + 470, y0 + 152), fill=BLUE, width=4, head=16)
    text(draw, (x1 - 212, y0 + 104), "70° cone half-angle", size=22, fill=BLUE, anchor="la")


def draw_assumptions(draw: ImageDraw.ImageDraw, box):
    draw.rounded_rectangle(box, radius=28, fill=LIGHT_BLUE, outline=(176, 205, 244), width=3)
    x0, y0, _, _ = box
    text(draw, (x0 + 28, y0 + 26), "Current aerodynamic assumptions", size=30, bold=True, fill=BLUE)
    text(
        draw,
        (x0 + 30, y0 + 82),
        "• Body flown at 0° nominal AoA\n"
        "• Panels shown in the symmetric broadside drag state\n"
        "• Differential left/right deployment appears only in exploratory crossrange studies",
        size=23,
        fill=TEXT,
        spacing=12,
    )


def main():
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)

    img = Image.new("RGB", (W, H), BG)
    draw = ImageDraw.Draw(img)
    draw.rounded_rectangle((48, 44, W - 48, H - 44), radius=34, outline=OUTLINE, width=3, fill=CARD)

    text(draw, (92, 88), "Hard-Lander States And Analysis Geometry", size=54, bold=True, fill=NAVY)
    text(
        draw,
        (95, 152),
        "Simplified communication figure for the Mars hard-lander study. This is an analysis-consistent surrogate, not exact SHIELD flight hardware.",
        size=24,
        fill=MUTED,
    )

    draw_state_card(
        draw,
        (92, 232, 1060, 712),
        title="Baseline body-only state",
        subtitle="Passive hard-lander configuration used as the high-β reference.",
        deployed=False,
    )
    draw_state_card(
        draw,
        (92, 748, 1060, 1228),
        title="Symmetric deployed drag state",
        subtitle="Both panels deployed equally to reduce ballistic coefficient.",
        deployed=True,
    )
    draw_top_view(draw, (1120, 248, 2106, 620))
    draw_side_geometry(draw, (1120, 652, 2106, 988))
    draw_assumptions(draw, (1120, 1012, 2106, 1228))

    img.save(OUT_PATH)
    VAULT_COPY.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(OUT_PATH, VAULT_COPY)
    print(OUT_PATH)


if __name__ == "__main__":
    main()
