# create_overlapping_tiles warns and errors on invalid inputs

    Code
      tiles <- create_overlapping_tiles(r, overlap = 200, sz = 10000)
    Condition
      Warning in `create_overlapping_tiles()`:
      Requested tile size larger than input area - returning a single tile of whole area!

---

    Code
      create_overlapping_tiles(r, overlap = 5000, sz = 5000)
    Condition
      Error in `create_overlapping_tiles()`:
      ! overlap must be zero or positive and less than sz!!

---

    Code
      tiles <- create_overlapping_tiles(r, overlap = 200, sz = 2025)
    Condition
      Warning in `create_overlapping_tiles()`:
      Choice of tile size is NOT divisible by resolution of template.r!!

