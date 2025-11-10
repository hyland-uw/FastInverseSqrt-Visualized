# Data Preparation

The scripts in this directory regenerate the data objects that ship with the
`visualfrsr` package. They read the canonical CSV exports under `inst/extdata/`
and convert them into compressed `.rda` files stored in `data/`.

To recreate the bundled datasets:

```sh
Rscript data-raw/prepare_binned_data.R
```

Run these scripts whenever the upstream generators in `inst/extdata/` change.
