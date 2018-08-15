# double data.table while inter-switching specific columns; e.g. when complementing DMS datatable such that Pos1 and Pos2 are symmetric
switch_double_DT = function(DT,cols_switchdouble,cols_double) {
  # cols_switchdouble: list of columns pairs
  sd_text = ""
  for (i in 1:length(cols_switchdouble)) {
    sd_text = c(sd_text,paste0(cols_switchdouble[[i]][1]),"=c(",cols_switchdouble[[i]][1],",",cols_switchdouble[[i]][2],"),")
    sd_text = c(sd_text,paste0(cols_switchdouble[[i]][2]),"=c(",cols_switchdouble[[i]][2],",",cols_switchdouble[[i]][1],"),")
  }
  eval(parse(text = paste0("DT = DT[,.(",paste0(sd_text,collapse=""),paste0(cols_double,collapse=","),")]")))
}
