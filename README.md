# OpenGDC
Python implementation of the OpenGDC software https://github.com/fabio-cumbo/OpenGDC

### Usage

```
python OpenGDC.py [--tumor          [GDC_TUMOR]             ]
                  [--datatype       [GDC_DATATYPE]          ]
                  [--after          [AFTER_DATETIME]        ]
                  [--download       [DOWNLOAD_FLAG]         ]
                  [--download_dir   [DOWNLOAD_DIRECTORY]    ]
                  [--convert        [CONVERT_FLAG]          ]
                  [--convert_dir    [CONVERT_DIRECTORY]     ]
                  [--matrix         [EXPORT_TO_MATRIX]      ]
                  [--settings       [SETTINGS_FILE]         ]
                  [--verbose        [VERBOSE_FLAG]          ]

optional arguments:
    --after  [AFTER_DATETIME]              default value None
    --matrix [EXPORT_TO_MATRIX]            default value 0
    --dimensionality [HD_DIMENSION]       default value 10000
    --retrain [RETRAINING_ITERATIONS]     default value 1

Both --tumor and --datatype are case sensitive;
Flags --download and --convert activate the download and conversion of the GDC data;
The convert procedure requires both the --download_dir and --convert-dir, additionally to the --convert flag enabled.
```

### Credits

Please credit our work in your manuscript by citing:

> Eleonora Cappelli, Fabio Cumbo, Anna Bernasconi, Arif Canakoglu, Stefano Ceri, Marco Masseroli, and Emanuel Weitschek. "OpenGDC: unifying, modeling, integrating cancer genomic data and clinical metadata" Appl. Sci. 2020, 10(18), 6367. [https://doi.org/10.3390/app10186367](https://doi.org/10.3390/app10186367)
