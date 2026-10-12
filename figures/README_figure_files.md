# Figure files — which ones to use

There are two sets of PNG renders in this directory. **Only the second set is suitable for
submission.**

## `figures/*.png` — screen-resolution previews, DO NOT SUBMIT

These are small (900–1300 px wide) and carry **no dpi metadata**. Printed at double-column
width (180 mm) they resolve to roughly 120–190 dpi and would fail the journal's ≥300 dpi
requirement. They are kept only for quick inspection on screen.

## `figures/SR_submission/*.png` — submission-quality renders, USE THESE

| Property | Value |
|---|---|
| Pixel width | 4323 px (all seven figures) |
| Embedded resolution | 600 dpi |
| Effective resolution at 180 mm width | 610 dpi |
| Font sizes | 10–17 pt in the SVG sources; smallest used is 10 pt, above the 6 pt minimum |

All seven figures (Fig1–Fig7) are present at this resolution.

## Source files

`figures/*.svg` are the vector sources and should be preferred where the journal accepts
vector artwork; they contain no rasterisation artefacts and scale to any size.
