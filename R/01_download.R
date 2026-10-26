source(here::here("R", "00_packages.R"))
source(here::here("R", "00_config.R"))

.in  <- here::here("R", "01_download.R")
.out <- here::here("data", "raw", "download_manifest.csv")
if (file.exists(.out) && file.mtime(.out) > file.mtime(.in)) {
  message("  Up to date, skipping")
} else {

tables <- c("DEMO", "DPQ", "HSCRP", "BMX", "GHB", "TCHOL", "HDL",
            "TRIGLY", "BPX", "BPXO", "CBC", "SMQ", "ALQ", "PAQ",
            "DIQ", "MCQ", "HIQ")

dir.create(here::here("data", "raw"), recursive = TRUE, showWarnings = FALSE)
message("  Downloading ", nrow(cycles), " NHANES cycles x ", length(tables), " tables...")

for (i in seq_len(nrow(cycles))) {
  cyc <- cycles[i, ]
  message("--- ", cyc$label, " ---")
  for (tbl in tables) {
    nm <- paste0(tbl, "_", cyc$id)
    f <- here::here("data", "raw", paste0(nm, ".rds"))

    if (file.exists(f)) {
      tryCatch({ readRDS(f); message("  ", nm, " exists") },
               error = function(e) { file.remove(f); message("  ", nm, " corrupt") })
      if (file.exists(f)) next
    }

    tryCatch({
      d <- nhanes(nm, translate = FALSE)
      if (!is.data.frame(d) || nrow(d) == 0) {
        message("  ", nm, " EMPTY (skipped)")
        next
      }
      saveRDS(d, f)
      message("  ", nm, " OK (", nrow(d), " rows)")
    }, error = function(e) {
      # P_ fallback for _J tables that may only exist under pre-pandemic naming
      if (cyc$id == "J") {
        p_nm <- paste0("P_", tbl)
        message("  ", nm, " failed, trying ", p_nm, "...")
        tryCatch({
          d <- nhanes(p_nm, translate = FALSE)
          if (is.data.frame(d) && nrow(d) > 0) {
            saveRDS(d, f)
            message("  ", p_nm, " -> ", nm, " OK (", nrow(d), " rows) [P_ fallback]")
          } else {
            message("  ", p_nm, " EMPTY (skipped)")
          }
        }, error = function(e2) message("  ", p_nm, " FAILED: ", e2$message))
      } else {
        message("  ", nm, " FAILED: ", e$message)
      }
    })
  }
}

# --- MD5 manifest for data drift detection ---
raw <- list.files(here::here("data", "raw"), pattern = "^[A-Z].*\\.rds$",
                  full.names = TRUE)
new_manifest <- data.frame(
  table    = tools::file_path_sans_ext(basename(raw)),
  md5      = vapply(raw, tools::md5sum, "", USE.NAMES = FALSE),
  snapshot = Sys.Date(),
  nhanesA  = as.character(packageVersion("nhanesA"))
)
if (file.exists(here::here("data", "raw", "download_manifest.csv"))) {
  old <- read.csv(here::here("data", "raw", "download_manifest.csv"),
                  stringsAsFactors = FALSE)
  shared <- intersect(old$table, new_manifest$table)
  changed <- shared[old$md5[match(shared, old$table)] !=
                    new_manifest$md5[match(shared, new_manifest$table)]]
  if (length(changed) > 0) {
    msg <- paste("Data files changed since last snapshot:",
                 paste(changed, collapse = ", "),
                 "\nSet ALLOW_DATA_DRIFT=TRUE to override.")
    if (identical(Sys.getenv("ALLOW_DATA_DRIFT"), "TRUE")) warning(msg) else stop(msg)
  }
}
write.csv(new_manifest, here::here("data", "raw", "download_manifest.csv"),
          row.names = FALSE)
message("  Manifest: ", nrow(new_manifest), " tables checksummed.")
}
