#collection of scripts used in DMS2struct project


# convert amino acid abbreviations between 1 and 3 letter code
convert_AAabr_one_three =  function(input,arg.len=length(unlist(strsplit(input,"")))) {
  one = unlist(strsplit("ACDEFGHIKLMNPQRSTVWY",""))
  a ="AlaCysAspGluPheGlyHisIleLysLeuMetAsnProGlnArgSerThrValTrpTyr"
  three = sapply(seq(1,nchar(a),by=3), function(x){substr(a, x, x+3-1)})
  if (arg.len == 1) {
    out = three[one %in% toupper(input)]
  } else if (arg.len == 3) {
    if (toupper(input) == "TER") {
      out = input
    } else {
      out = one[toupper(three) %in% toupper(input)]
    }
  } else {
    out = "XXX"
  }
  return(out)
}



