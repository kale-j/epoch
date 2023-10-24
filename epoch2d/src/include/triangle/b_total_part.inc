        bx_part = &
              hy(-1) * (gx(-1) * bx_total(cell_x1-1,cell_y2-1) &
            +           gx( 0) * bx_total(cell_x1  ,cell_y2-1) &
            +           gx( 1) * bx_total(cell_x1+1,cell_y2-1)) &
            + hy( 0) * (gx(-1) * bx_total(cell_x1-1,cell_y2  ) &
            +           gx( 0) * bx_total(cell_x1  ,cell_y2  ) &
            +           gx( 1) * bx_total(cell_x1+1,cell_y2  )) &
            + hy( 1) * (gx(-1) * bx_total(cell_x1-1,cell_y2+1) &
            +           gx( 0) * bx_total(cell_x1  ,cell_y2+1) &
            +           gx( 1) * bx_total(cell_x1+1,cell_y2+1))

        by_part = &
              gy(-1) * (hx(-1) * by_total(cell_x2-1,cell_y1-1) &
            +           hx( 0) * by_total(cell_x2  ,cell_y1-1) &
            +           hx( 1) * by_total(cell_x2+1,cell_y1-1)) &
            + gy( 0) * (hx(-1) * by_total(cell_x2-1,cell_y1  ) &
            +           hx( 0) * by_total(cell_x2  ,cell_y1  ) &
            +           hx( 1) * by_total(cell_x2+1,cell_y1  )) &
            + gy( 1) * (hx(-1) * by_total(cell_x2-1,cell_y1+1) &
            +           hx( 0) * by_total(cell_x2  ,cell_y1+1) &
            +           hx( 1) * by_total(cell_x2+1,cell_y1+1))

        bz_part = &
              hy(-1) * (hx(-1) * bz_total(cell_x2-1,cell_y2-1) &
            +           hx( 0) * bz_total(cell_x2  ,cell_y2-1) &
            +           hx( 1) * bz_total(cell_x2+1,cell_y2-1)) &
            + hy( 0) * (hx(-1) * bz_total(cell_x2-1,cell_y2  ) &
            +           hx( 0) * bz_total(cell_x2  ,cell_y2  ) &
            +           hx( 1) * bz_total(cell_x2+1,cell_y2  )) &
            + hy( 1) * (hx(-1) * bz_total(cell_x2-1,cell_y2+1) &
            +           hx( 0) * bz_total(cell_x2  ,cell_y2+1) &
            +           hx( 1) * bz_total(cell_x2+1,cell_y2+1))
