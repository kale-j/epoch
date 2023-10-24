        bx_part = &
              gx(-1) * bx_total(cell_x1-1) &
            + gx( 0) * bx_total(cell_x1  ) &
            + gx( 1) * bx_total(cell_x1+1)

        by_part = &
              hx(-1) * by_total(cell_x2-1) &
            + hx( 0) * by_total(cell_x2  ) &
            + hx( 1) * by_total(cell_x2+1)

        bz_part = &
              hx(-1) * bz_total(cell_x2-1) &
            + hx( 0) * bz_total(cell_x2  ) &
            + hx( 1) * bz_total(cell_x2+1)
