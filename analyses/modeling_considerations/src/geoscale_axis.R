library(geoscale)

geoscaleAxis = function(timescale,
                        side     = 1,
                        fraction = 0.8,
                        units    = c("Epoch", "Period"),
                        las.unit = c(2, 1),
                        cex.lab  = 0.7,
                        lwd.axis = 1,
                        las.axis = 1,
                        cex.axis = 0.7,
                        outer    = FALSE) {

  # compute the unit fraction (the amount of the margin to take up)
  units_end   = cumsum(las.unit) / sum(las.unit) * fraction
  units_start = c(0, units_end[-length(units_end)])

  # plot each unit
  num_units = length(units)
  for(i in 1:num_units) {

    # get this unit
    this_unit_type = units[i]
    this_unit = timescale[timescale$Type == this_unit_type,]

    # get the coordinates (in fractions for the marginal dimension)
    if ( side %in% c(1,3) ) {
      x0 = this_unit$Start
      x1 = this_unit$End
      y0 = units_start[i]
      y1 = units_end[i]
    } else if ( side %in% c(2,4) ) {
      stop("Not yet implemented.")
    } else {
      stop("Invalid side. Choose 1, 2, 3 or 4.")
    }

    # discard any unplottable units
    if ( side %in% c(1,3) ) {

      # determine the limit
      lim = par()$usr[1:2]

      # remove any units that end before the upper limit
      this_unit = this_unit[this_unit$End < max(lim),]

      # truncate any start points to the upper limit
      this_unit$Start = pmin(this_unit$Start, max(lim))

      # recompute the x limits
      x0 = this_unit$Start
      x1 = this_unit$End

    } else if ( side %in% c(2,4) ) {
      stop("Not yet implemented.")
    } else {
      stop("Invalid side. Choose 1, 2, 3 or 4.")
    }

    # compute the actual coordinates
    if ( side == 1 ) {
        if ( outer == FALSE ) {
            y0 = par()$mai[1] - y0 * par()$mai[1]
            y1 = par()$mai[1] - y1 * par()$mai[1]
        } else {
            y0 = par()$omi[1] - y0 * par()$omi[1]
            y1 = par()$omi[1] - y1 * par()$omi[1]
        }
        y0 = grconvertY(y0, from="inches", to="user")
        y1 = grconvertY(y1, from="inches", to="user")
    } else if ( side == 2 ) {
      stop("Not yet implemented.")
    } else if ( side == 3 ) {
      y0 = par()$mai[3] + y0 * par()$mai[3]
      y1 = par()$mai[3] + y1 * par()$mai[3]
      y0 = grconvertY(y0, from="inches", to="user")
      y1 = grconvertY(y1, from="inches", to="user")
    } else if ( side == 4 ) {
      stop("Not yet implemented.")
    } else {
      stop("Invalid side. Choose 1, 2, 3 or 4.")
    }

    # make the colors per unit
    unit_colors = rgb(this_unit$Col_R, this_unit$Col_G, this_unit$Col_B, maxColorValue = 255)

    # check the size of unit names
    unit_names = as.character(this_unit$Name)
    if ( las.unit[i] == 2 ) {

      # get the dimensions
      str_heights = strwidth(unit_names, units="inches", cex=cex.lab) * 1.2
      str_widths  = strheight(unit_names, units="inches", cex=cex.lab) * 1.2

      # check if they are within the box
      dont_fit = str_heights > abs(grconvertY(y1, "user", "inches") - grconvertY(y0, "user", "inches"))
      unit_names[dont_fit] = as.character(this_unit$Abbrev[dont_fit])
      str_heights = strwidth(unit_names, units="inches", cex=cex.lab) * 1.2

      # determine what to plot
      include_x = str_widths  < abs(grconvertX(x1, "user", "inches") - grconvertX(x0, "user", "inches"))
      include_y = str_heights < abs(grconvertY(y1, "user", "inches") - grconvertY(y0, "user", "inches"))

      # record the angle
      srt = 90

    } else if ( las.unit[i] == 1 ) {

      # get the dimensions
      str_heights = strheight(unit_names, units="inches", cex=cex.lab) * 1.2
      str_widths  = strwidth(unit_names, units="inches", cex=cex.lab) * 1.2

      # check if they are within the box
      dont_fit = str_widths > abs(grconvertX(x1, "user", "inches") - grconvertX(x0, "user", "inches"))
      unit_names[dont_fit] = as.character(this_unit$Abbrev[dont_fit])
      str_widths = strheight(unit_names, units="inches", cex=cex.lab) * 1.2

      # determine what to plot
      include_x = str_widths  < abs(grconvertX(x1, "user", "inches") - grconvertX(x0, "user", "inches"))
      include_y = str_heights < abs(grconvertY(y1, "user", "inches") - grconvertY(y0, "user", "inches"))

      srt = 0

    } else {
      stop("Invalid las. Choose 1 or 2.")
    }

    # determine coordinates
    ymid = grconvertY((grconvertY(y0, "user", "inches") + grconvertY(y1, "user", "inches")) / 2, "inches", "user")
    xmid = grconvertX((grconvertX(x0, "user", "inches") + grconvertX(x1, "user", "inches")) / 2, "inches", "user")
    if ( length(xmid) == 1 ) {
      xmid = rep(xmid, nrow(this_unit))
    }
    if ( length(ymid) == 1 ) {
      ymid = rep(ymid, nrow(this_unit))
    }

    # plot units
    rect(xleft=x0, xright=x1, ybottom=y0, ytop=y1, xpd=NA, col=unit_colors)

    # plot unit names (that fit within boxes)
    text(x=xmid[include_x & include_y], y=ymid[include_x & include_y], labels=unit_names[include_x & include_y], xpd=NA, srt=srt, cex=cex.lab)

  }

  # add axis at the bottom
  if ( side == 1 ) {
      if ( outer == FALSE ) {
          axis_line = (par()$mai[side] - grconvertY(y1, "user", "inches")) / par('cin')[2] * par('cex') * par('lheight')
      } else {
          axis_line = par()$oma[side] * fraction
      }
      axis(at=this_unit$Start[-length(this_unit$Start)], side=side, line=axis_line, lwd=0, lwd.tick=lwd.axis, las=las.axis, cex.axis=cex.axis)
  }

}
