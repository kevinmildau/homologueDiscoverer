#' Title
#'
#' @param peaks
#' @param fileTo
#' @param cid
#' @param crt
#' @param cmz
#' @param crtmin
#' @param crtmax
#' @param cmzmin
#' @param cmzmax
#' @param cquality
#'
#' @return
#' @export
#'
#' @examples
exportAsFeatureXML <- function(peaks, fileTo, cid = "Num", crt = "rt",
                               cmz = "mz", crtmin = "rtmin", crtmax = "rtmax",
                               cmzmin = "mzmin", cmzmax = "mzmax",
                               cquality = "sn"){

  lines = c()
  lines = c(lines, sprintf('<?xml version="1.0" encoding="ISO-8859-1"?>'))
  lines = c(lines, sprintf('  <featureMap version="1.4" id="fm_16311276685788915066" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/FeatureXML_1_4.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'))
  lines = c(lines, sprintf('    <dataProcessing completion_time="2011-02-09T11:23:14">'))
  lines = c(lines, sprintf('      <software name="XCMS" version="3.4.4" />'))
  lines = c(lines, sprintf('    </dataProcessing>'))
  lines = c(lines, sprintf('    <featureList count="%d">', nrow(peaks)))
  fileConn<-file(fileTo, open="w")
  writeLines(lines, fileConn)
  close(fileConn)

  for(rowi in 1:nrow(peaks)){

    id = ifelse(cid %in% colnames(peaks), peaks[rowi, cid], sprintf("Row%s", rownames(peaks)[rowi]))
    rt = peaks[rowi, crt]
    mz = peaks[rowi, cmz]
    rtmin = ifelse(crtmin %in% colnames(peaks), peaks[rowi, crtmin], rt - 5.)
    rtmax = ifelse(crtmax %in% colnames(peaks), peaks[rowi, crtmax], rt + 5.)
    mzmin = ifelse(cmzmin %in% colnames(peaks), peaks[rowi, cmzmin], mz * (1-5/1E6))
    mzmax = ifelse(cmzmax %in% colnames(peaks), peaks[rowi, cmzmax], mz * (1+5/1E6))
    quality = ifelse(cquality %in% colnames(peaks), peaks[rowi, cquality], -1)

    lines = c()
    lines = c(lines, sprintf('<feature id="%s">', id))
    lines = c(lines, sprintf('  <position dim="0">%f</position>', rt))
    lines = c(lines, sprintf('  <position dim="1">%f</position>', mz))
    lines = c(lines, sprintf('  <intensity>%f</intensity>', 1.))
    lines = c(lines, sprintf('  <quality dim="0">0</quality>'))
    lines = c(lines, sprintf('  <quality dim="1">0</quality>'))
    lines = c(lines, sprintf('  <overallquality>%s</overallquality>', quality))
    lines = c(lines, sprintf('  <charge>1</charge>'))
    lines = c(lines, sprintf('  <convexhull nr="0">'))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmin, mzmin))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmin, mzmax))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmax, mzmin))
    lines = c(lines, sprintf('    <pt x="%f" y="%f" />', rtmax, mzmax))
    lines = c(lines, sprintf('  </convexhull>'))
    lines = c(lines, sprintf('</feature>'))
    fileConn<-file(fileTo, open="a")
    writeLines(lines, fileConn)
    close(fileConn)
  }

  lines = c()
  lines = c(lines, sprintf('    </featureList>'))
  lines = c(lines, sprintf('  </featureMap>'))
  fileConn<-file(fileTo, open="a")
  writeLines(lines, fileConn)
  close(fileConn)
}
