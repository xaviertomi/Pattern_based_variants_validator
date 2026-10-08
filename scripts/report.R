args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop('Usage: Rscript scripts/report.R OUTPUT_DIR')
root <- normalizePath(args[[1]], mustWork = TRUE)
reports <- file.path(root, 'reports')
dir.create(reports, recursive = TRUE, showWarnings = FALSE)
read_table <- function(name) {
  path <- file.path(reports, paste0(name, '.tsv'))
  if (!file.exists(path)) stop(paste('Missing report input:', path))
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
}
page <- function(title, text) {
  plot.new(); title(main = title)
  text(0.05, 0.85, paste(strwrap(text, width = 90), collapse = '\n'), adj = c(0, 1), cex = 0.9)
}
coverage <- read_table('motif_coverage')
decisions <- read_table('selection_decisions')
representatives <- read_table('representative_validation')
annotations <- read_table('annotation_counts')
number <- function(table, choices) {
  name <- choices[choices %in% names(table)]
  if (!length(name)) return(rep(NA_real_, nrow(table)))
  suppressWarnings(as.numeric(table[[name[[1]]]]))
}
pdf(file.path(reports, 'motif_selection.pdf'), width = 10, height = 7)
page('Sequence support, not variant truth', 'Family stability summarizes related models across all five held-out rotations. Exact representative validation is a separate matrix-specific test on its source rotation held-out fold across the predetermined required sample groups. External motif matches are sequence support, not variant truth; no match is unclassified, not false. FIMO p-values and Tomtom q-values are not probabilities that a variant is true. Overlapping nearby windows are non-independent; these summaries are descriptive and have no independent-window confidence intervals.')
ratio <- number(coverage, c('coverage_ratio', 'ratio'))
central <- number(coverage, c('central_proportion_difference', 'central_difference'))
real <- number(coverage, c('real_coverage', 'coverage_real'))
valid <- !is.na(ratio) & !is.na(central) & !is.na(real)
if (any(valid)) {
  finite <- ratio[valid & is.finite(ratio)]
  ceiling <- if (length(finite)) max(1, finite) * 1.1 else 1
  x <- ratio; x[is.infinite(x)] <- ceiling
  plot(x[valid], central[valid], pch = 21, bg = '#3175aa', cex = 3 * sqrt(pmax(0, real[valid])), xlab = 'Real/control window coverage ratio (Inf in capped band)', ylab = 'Conditional central proportion: real minus control', main = 'Held-out enrichment and conditional positions')
  abline(h = 0, v = 1, lty = 2, col = 'gray')
  if (any(valid & is.infinite(ratio))) text(x[valid & is.infinite(ratio)], central[valid & is.infinite(ratio)], labels = 'Inf', pos = 4, cex = 0.7)
  mtext('Point AREA is proportional to real coverage; TSV values remain unchanged.', side = 3, line = 0.2, cex = 0.8)
} else page('No defined enrichment/central comparison', 'No rows have both a defined enrichment ratio and a defined conditional central difference. Zero-hit conditional proportions remain undefined, never zero. An empty shortlist is a valid scientific outcome.')
page('Separate eligibility gates', sprintf('Recorded family decisions: %d. Recorded exact-representative dataset evaluations: %d. Every retained family must pass family stability AND have an original exact matrix passing the cross-group representative gate. A family union cannot substitute for a representative matrix. Thresholds are configurable sequence-support heuristics, not calibrated variant-truth criteria.', nrow(decisions), nrow(representatives)))
if (nrow(annotations)) {
  totals <- tapply(annotations$record_count, annotations$annotation_status, sum)
  barplot(totals, las = 2, col = '#3175aa', ylab = 'Original genotype rows', main = 'External annotation statuses')
} else page('No external genotype rows', 'The external input yielded no genotype rows. Output tables retain their complete column schema.')
dev.off()
pdf(file.path(reports, 'conditional_positions.pdf'), width = 10, height = 7)
position_path <- file.path(reports, 'conditional_position_bins.tsv')
if (file.exists(position_path)) {
  positions <- read.delim(position_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c('', 'NA'))
  bins <- number(positions, c('midpoint_bin', 'bin'))
  difference <- number(positions, c('conditional_difference', 'real_minus_control', 'difference'))
  good <- is.finite(bins) & is.finite(difference)
  if (any(good)) {
    plot(bins[good], difference[good], pch = 16, xlab = 'floor(reference-oriented midpoint distance), bp', ylab = 'Conditional hit-window proportion: real minus control', main = 'Highest-score occurrence positions')
    abline(h = 0, v = 0, lty = 2, col = 'gray')
  } else page('Undefined conditional positions', 'No defined real-minus-control conditional position proportions are available. Zero-hit controls have undefined conditional distributions and are never replaced by zero.')
} else page('No selected conditional position comparison', 'No selected motif with defined real and control conditional positional evidence is available. Positions use reference orientation on both strands and highest-score occurrences; one-base bins are floor(midpoint distance), including half-base centers.')
page('Position interpretation', 'Conditional proportions use only motif-containing windows, not all genotype records. The highest numeric-score occurrence represents each motif/window, with reported p-value and deterministic coordinate tie breaks; it is not the closest-to-center occurrence. Exact representative evidence is source-rotation-specific, not five independent representative validations.')
dev.off()
writeLines(c('# Motif workflow report', '', 'This report describes sequence support only. A promising external ALT genotype has a reference-context motif match; it is not established as a true variant. Unclassified means no qualifying sequence match, never false.', '', '## Two separate selection gates', 'Family stability measures related, transferable members over five whole-contig held-out rotations. Exact representative validation evaluates each original candidate matrix on its source rotation held-out fold in every predetermined required PacBio group and every mapped Illumina dataset. Only matrices passing the second gate can freeze. The representative is neither refitted nor a family union.', '', '## Descriptive evidence', 'FIMO reported p-values and Tomtom q-values are not variant-truth probabilities. Nearby windows overlap, shared reference contexts across technologies are non-independent, and positional overlap is not normalized allele concordance. No independent-window tests or confidence intervals are claimed.', 'Point area in motif_selection.pdf is proportional to real window coverage. Infinite enrichment occupies a labelled capped display band only; TSV values retain Inf. Conditional position distributions include only hit windows, use the highest-score occurrence in reference orientation and floor(midpoint distance) one-base bins. Undefined zero-hit control distributions stay undefined.', '', '## External identity and annotations', 'Inspect annotation_counts.tsv and application/manifest.json for identity_check_status. header_names_only explicitly means biological independence and aliases have not been verified. All original genotypes remain in application/genotypes.tsv.gz; application/alt_annotations.tsv.gz contains ALT carriers and all_matches.tsv.gz retains every qualifying occurrence.', 'Reference-only and entirely missing calls are not_applicable; extraction failures are not_scanned; matches are promising; scan-eligible no-hit calls are unclassified. Empty libraries produce unclassified/no_selected_motifs for eligible ALT contexts without inventing scan failures.', '', '## Evidence tables', 'dataset_fold_counts.tsv, exclusions.tsv, streme_usage.tsv, motif_coverage.tsv, transfer_evidence.tsv, family_coverage.tsv, family_membership.tsv, selection_decisions.tsv, representative_validation.tsv and annotation_counts.tsv retain the denominators, failures and separate gate decisions. Scientific zero-hit and empty outcomes are valid; tool failures abort rather than becoming zero support.'), file.path(reports, 'README.md'), useBytes = TRUE)
