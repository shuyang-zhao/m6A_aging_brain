# Tints and shades of a colour: a tint fills an area, the colour itself strokes it, so that a hue reads the same at
# every mark size
tint <- function(col, amount = 0.6) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(m + (1 - m) * amount))
}
shade <- function(col, amount = 0.3) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(m * (1 - amount)))
}
