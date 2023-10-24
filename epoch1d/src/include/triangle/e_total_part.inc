        ex_part = &
              hx(-1) * ex_total(cell_x2-1) &
            + hx( 0) * ex_total(cell_x2  ) &
            + hx( 1) * ex_total(cell_x2+1)

        ey_part = &
              gx(-1) * ey_total(cell_x1-1) &
            + gx( 0) * ey_total(cell_x1  ) &
            + gx( 1) * ey_total(cell_x1+1)

        ez_part = &
              gx(-1) * ez_total(cell_x1-1) &
            + gx( 0) * ez_total(cell_x1  ) &
            + gx( 1) * ez_total(cell_x1+1)
