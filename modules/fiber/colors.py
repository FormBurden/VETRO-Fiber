# modules/fiber/colors.py

# Canonical 1–12 color order (rows of 12), per project rules:
# 1-Blue, 2-Orange, 3-Green, 4-Brown, 5-Slate, 6-White,
# 7-Red, 8-Black, 9-Yellow, 10-Violet, 11-Rose, 12-Aqua
FIBER_COLORS = [
    "Blue","Orange","Green","Brown","Slate","White",
    "Red","Black","Yellow","Violet","Rose","Aqua"
]

def fiber_num_to_color_label(fiber_num: int) -> str:
    """
    Map an absolute fiber number (1..N) to the canonical 1..12 color index + name.
    Examples:
      12 -> "12 - Aqua"
      13 -> "1 - Blue"
      15 -> "3 - Green"
    """
    idx = ((int(fiber_num) - 1) % 12) + 1
    return f"{idx} - {FIBER_COLORS[idx-1]}"
