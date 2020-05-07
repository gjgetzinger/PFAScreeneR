#' Make a PFAS structure database from a list of input files
#'
#' Reads input files, performs transformations, writes results to database
#' (optionally keeps intermediate sdf files)
#'
#' @param mollist_path A path to a directory containing .xls molecule lists for
#'   import
#' @param biotrans_jarfile A path to a BioTransformer jar file
#' @param result_path A directory to store results (defaults to currect working
#'   directory)
#' @param keep_sdf T/F should intermediate sdf files be preserved
#'
#' @return Filenames were results are written
#' @export
#'
make_pfas_db <- function(
  mollist_path = NULL,
  biotrans_jarfile = NULL,
  result_path = NULL,
  keep_sdf = F
){
  stopifnot(reticulate::py_available(initialize = T))
  py_mods <- c(
    "rdkit","pandas","numpy",
    "molvs","sys","os","subprocess",
    "tempfile","fnmatch","xlrd",
    "sqlite3","progress","time"
    )
  stopifnot(sapply(py_mods, reticulate::py_module_available))

  py_file <- system.file(package = 'PFAScreeneR') %>%
    list.files(full.names = T,recursive = T,
               pattern = 'pfas_masslist.py')

  py_file <- normalizePath(py_file, mustWork = T)

  reticulate::source_python(file = py_file)

  biotrans_jarfile <- normalizePath(biotrans_jarfile, mustWork = T)
  biotrans_required_files <-
    list.dirs(dirname(biotrans_jarfile),
              recursive = F,
              full.names = F)

  if (!all(c('database', 'supportfiles') %in%  biotrans_required_files)) {
    stop(
      "BioTransformer requires database and supportfiles in
    the same directory as the jar file.
    See: https://bitbucket.org/djoumbou/biotransformerjar/src/master/
    for installation instructions"
    )
  }

  result_path <-
    paste0(normalizePath(result_path, mustWork = T), '/PFAScreenR%s')

  result_db <- sprintf(result_path, '.db')
  hyd_sdf <- sprintf(result_path, '_hyd_out.sdf')
  biotrans_in_sdf <- sprintf(result_path, '_biotrans_in.sdf')
  biotrans_out_sdf <- sprintf(result_path, '_biotrans_out.sdf')
  result_sdf <- sprintf(result_path, '_result.sdf')

  make_pfas_masslist(
    mollist_path = mollist_path,
    biotrans_jarfile = biotrans_jarfile,
    hyd_sdf = hyd_sdf,
    biotrans_in_sdf = biotrans_in_sdf,
    biotrans_out_sdf = biotrans_out_sdf,
    result_sdf = result_sdf,
    result_db = result_db
  )

  if (!keep_sdf) {
    unlink(x = c(hyd_sdf, biotrans_in_sdf, biotrans_out_sdf, result_sdf))
    return(result_db)
  } else {
    return(c(
      result_db,
      hyd_sdf,
      biotrans_in_sdf,
      biotrans_out_sdf,
      result_sdf
    ))
  }
}

#' Make a flat file from the PFAScreeneR.db file
#'
#' Useful for creating flat files to use as input to other software
#' (e.g., To make a mass list in Thermo Compound Discoverer)
#' @param db_file Full path to a database file
#' @param tables List of table names to merge and write to file
#'
#' @return filename of the exported flat file
#' @export
#'
make_pfas_flat <-
  function(db_file,
           tables = c('molecules',
                      'list_membership',
                      'mol_props')) {
    db_file <- normalizePath(db_file, mustWork = T)
    conn <- DBI::dbConnect(RSQLite::SQLite(), db_file)
    tables <- match.arg(tables, several.ok = T)
    stopifnot(tables %in% DBI::dbListTables(conn))
    dat <- purrr::map(tables, .f = function(x) dplyr::tbl(src = conn, x))
    stopifnot(any(
      sapply(tables, function(x)
        'ID' %in% DBI::dbListFields(conn, x), USE.NAMES = F)
    ))

    x <- dat[[1]]
    for (y in dat[2:length(tables)]) {
      x <- dplyr::left_join(x, y)
    }

    dplyr::as_tibble(x) %>%
    dplyr::group_by(ID) %>%
      dplyr::summarise_all(.funs = list(~paste(unique(.), collapse = ','))) %>%
      readr::write_csv(x = ., path = gsub('[.]db', '_dump.csv', db_file),col_names = T)
  }

globalVariables(c(".", "ID", "make_pfas_masslist"))
