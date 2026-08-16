* 2026.07.13
    * Take a look at implicit/explicit integration
    * For calibration
        * For tension, use direct tension test, without rebars
        * For compression, use direct compression, without rebars
        * Use 0.1m * 0.1m * 0.1 cube for calibration, beam-column model will do the work, both pure tension and compression
            * Use 1.5 horizon lattice for the same model
            * With tension part, find tensile parameters of the material model
            * Do something similar for compression
            * Inform Beyazit Hoca
        * For numerical part of the pushover, modified newton initial
