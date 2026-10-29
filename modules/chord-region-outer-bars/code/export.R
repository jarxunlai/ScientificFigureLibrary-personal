# 设备与通用几何工具；分组、范围、配色、数据选择保留在各 organized.R。
# 当前锁定环境缺少 svglite/rsvg，使用现有 ragg + Cairo；不自动装包。
options(stringsAsFactors = FALSE)
options(device = function(...) grDevices::pdf(file = NULL, ...))
save_base <- function(expr, stem, width_mm = 180, height_mm = 180, dpi = 300) {
  code <- substitute(expr)
  env <- parent.frame()
  for (format in c("png", "pdf")) {
    if (format == "png") {
      ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm,
                    units = "mm", res = dpi, background = "white")
    } else {
      grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm/25.4,
                          height = height_mm/25.4, family = "Arial", bg = "white")
    }
    on.exit(grDevices::dev.off(), add = TRUE)
    par(family = "Arial", ps = 7, mar = rep(0, 4), xpd = NA)
    eval(code, envir = env)
    grDevices::dev.off()
    on.exit(NULL, add = FALSE)
  }
}
canvas <- function(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3)) {
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, asp = 1, xaxs = "i", yaxs = "i")
}
xy <- function(angle, radius) {
  cbind(radius * cos(angle*pi/180), radius * sin(angle*pi/180))
}
arc <- function(start, end, radius, col = "grey40", lwd = .5, lty = 1) {
  a <- seq(start, end, length.out = max(8, ceiling(abs(end-start)*2)))
  lines(xy(a, radius), col = col, lwd = lwd, lty = lty)
}
sector <- function(start, end, inner, outer, col, border = NA, lwd = .4) {
  a <- seq(start, end, length.out = max(5, ceiling(abs(end-start)*2)))
  z <- rbind(xy(a, outer), xy(rev(a), inner))
  polygon(z, col = col, border = border, lwd = lwd)
}
radial_segment <- function(angle, from, to, col = "black", lwd = .5) {
  z <- xy(angle, c(from, to))
  segments(z[1,1], z[1,2], z[2,1], z[2,2], col = col, lwd = lwd)
}
radial_text <- function(angle, radius, label, col = "black", cex = 5/7,
                        facing = "radial", font = 1, outward = TRUE) {
  a <- (angle %% 360)
  z <- xy(angle, radius)
  left <- a > 90 & a < 270
  rotation <- if (facing == "tangent") a-90 else a
  if (facing == "tangent" && a > 180) rotation <- rotation + 180
  if (facing != "tangent" && left) rotation <- rotation + 180
  adj <- if (facing == "tangent") .5 else if (left == outward) 1 else 0
  text(z[1,1], z[1,2], label, srt = rotation, adj = c(adj,.5),
       col = col, cex = cex, font = font)
}
arc_text <- function(angle, radius, label, cex = 5/7, col = "black", font = 1) {
  chars <- strsplit(label, "", fixed = TRUE)[[1]]
  widths <- strwidth(chars, cex = cex, font = font, units = "user")
  widths[chars == " "] <- strwidth("n", cex = cex, units = "user")*.40
  angular_width <- widths/radius*180/pi
  direction <- if ((angle %% 360) > 180) 1 else -1
  positions <- angle + direction*(cumsum(angular_width)-angular_width/2-sum(angular_width)/2)
  for (i in seq_along(chars)) {
    p <- xy(positions[i],radius)
    rot <- positions[i]-90 + if(direction==1)180 else 0
    text(p[1,1],p[1,2],chars[i],cex=cex,col=col,font=font,srt=rot)
  }
}
curve_link <- function(a1, r1, a2, r2, control = .08, col = "grey70", lwd = .5) {
  p1 <- xy(a1,r1); p2 <- xy(a2,r2)
  c1 <- xy(a1,control); c2 <- xy(a2,control)
  tt <- seq(0,1,length.out = 80)
  bx <- (1-tt)^3*p1[1,1] + 3*(1-tt)^2*tt*c1[1,1] + 3*(1-tt)*tt^2*c2[1,1] + tt^3*p2[1,1]
  by <- (1-tt)^3*p1[1,2] + 3*(1-tt)^2*tt*c1[1,2] + 3*(1-tt)*tt^2*c2[1,2] + tt^3*p2[1,2]
  lines(bx, by, col = col, lwd = lwd)
}
raster_points <- function(x, y = NULL, pch = 16, col = "black", cex = .4) {
  coordinates <- xy.coords(x,y)
  usr <- par("usr"); inches <- par("pin")
  capture <- ragg::agg_capture(width=inches[1],height=inches[2],units="in",res=300,
                               background="transparent")
  par(mar=rep(0,4),ps=7,xaxs="i",yaxs="i")
  plot.new(); plot.window(xlim=usr[1:2],ylim=usr[3:4],xaxs="i",yaxs="i")
  points(coordinates$x,coordinates$y,pch=pch,col=col,cex=cex)
  layer <- capture(); dev.off()
  rasterImage(layer,usr[1],usr[3],usr[2],usr[4],interpolate=FALSE)
}
color_key <- function(fun, limits, x, y, width = .035, height = .35,
                       ticks = pretty(limits,3), title = "") {
  vals <- seq(limits[1],limits[2],length.out = 100)
  for (i in 1:99) rect(x,y+(i-1)*height/99,x+width,y+i*height/99,
                      col = fun(mean(vals[c(i,i+1)])), border = NA)
  text(x+width+.025,y+(ticks-limits[1])/diff(limits)*height,
       labels = format(ticks,trim = TRUE),adj = 0,cex = 5/7)
  text(x-.035,y+height/2,title,srt = 90,cex = 6/7)
}
