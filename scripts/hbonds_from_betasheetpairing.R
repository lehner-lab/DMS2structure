################################################################################
################ from beta sheet pairing derive hbond combinations #############
################################################################################
#subfunction to predict_beta_sheets function
hbonds_from_betasheetpairing = function(beta_sheet_pairing,ss_data = c()) {
  
  beta_hbonds = data.table(hn_opt1 = integer(),o_opt1 = integer(),hn_opt2 = integer(),o_opt2 = integer(),sheet= integer(),hbond = integer())
  # beta_sheet_pairing = copy(beta_sheet_pairing)
  if (nrow(beta_sheet_pairing) > 0) {
    for (i in 1:nrow(beta_sheet_pairing)) {
      L = beta_sheet_pairing[i,pos1_max - pos1_min + 1]
      if (beta_sheet_pairing[i,type=="anti-par"]) { #anti-parallel sheet
        
        if (L %% 2 == 1) { #if uneven, decide which side to keep; 
          
          if (length(ss_data) > 0) { #use the one wiht lower p_value
            p_start = ss_data[Pos1 == beta_sheet_pairing[i,pos1_min] & Pos2 == beta_sheet_pairing[i,pos2_max],beta_antipar_p]
            p_end = ss_data[Pos1 == beta_sheet_pairing[i,pos1_max] & Pos2 == beta_sheet_pairing[i,pos2_min],beta_antipar_p]
            if (p_start < p_end) {
              beta_sheet_pairing[i,':=' (pos1_max = pos1_max - 1, pos2_min = pos2_min + 1)]
            } else {
              beta_sheet_pairing[i,':=' (pos2_max = pos2_max - 1, pos1_min = pos1_min + 1)]
            }
          } else { #without additional info, keep lower position
            beta_sheet_pairing[i,':=' (pos2_max = pos2_max - 1, pos1_min = pos1_min + 1)]
          }
          L=L-1
        }
        beta_hbonds = rbind(beta_hbonds,beta_sheet_pairing[i,.(hn_opt1 = c(seq(pos1_min,pos1_max,2),seq(pos2_max,pos2_min,-2)),
                                                               o_opt1 = c(seq(pos2_max,pos2_min,-2),seq(pos1_min,pos1_max,2)),
                                                               hn_opt2 = c(seq(pos1_min+1,pos1_max,2),seq(pos2_max-1,pos2_min,-2)),
                                                               o_opt2 = c(seq(pos2_max-1,pos2_min,-2),seq(pos1_min+1,pos1_max,2)),
                                                               sheet=i,hbond=nrow(beta_hbonds)+1:L)])
        
      } else { #parallel sheet
        beta_hbonds = rbind(beta_hbonds,beta_sheet_pairing[i,.(hn_opt1 = c(seq(pos2_min+1,pos2_max,2),seq(pos1_min+2,pos1_max,2)),
                                                               o_opt1 = c(seq(pos1_min,pos1_max-1,2),seq(pos2_min+1,pos2_max-1,2)),
                                                               hn_opt2 = c(seq(pos1_min+1,pos1_max,2),seq(pos2_min+2,pos2_max,2)),
                                                               o_opt2 = c(seq(pos2_min,pos2_max-1,2),seq(pos1_min+1,pos1_max-1,2)),
                                                               sheet=i,hbond=nrow(beta_hbonds)+1:(L-1))])
      }
    }
    if (nrow(beta_sheet_pairing) > 1) {
      #create all possible combinations of sheets
      require(combinat)
      sheet_comb = list()
      for (i in 1:(nrow(beta_sheet_pairing)-1)) {
        x=combn(1:nrow(beta_sheet_pairing),i)
        sheet_comb = c(sheet_comb,split(x, rep(1:ncol(x), each = nrow(x))))
      }
      #swap options until hbonding is consistent
      iterations = 0
      while ((beta_hbonds[,.N,hn_opt1][,sum(N>1)>0] | beta_hbonds[,.N,hn_opt2][,sum(N>1)>0]) & iterations < 3) {
        for (i in 1:length(sheet_comb)) {
          beta_hbonds[sheet %in% sheet_comb[[i]], c("hn_opt1","o_opt1","hn_opt2","o_opt2") := .(hn_opt2,o_opt2,hn_opt1,o_opt1)]
          if (!beta_hbonds[,.N,hn_opt1][,sum(N>1)>0] & !beta_hbonds[,.N,hn_opt2][,sum(N>1)>0]) {
            break
          } else { #reverse
            beta_hbonds[sheet %in% sheet_comb[[i]], c("hn_opt1","o_opt1","hn_opt2","o_opt2") := .(hn_opt2,o_opt2,hn_opt1,o_opt1)]
          }
          if (i == length(sheet_comb)) {iterations = iterations + 1}
        }
      }
    }
  }
  return(beta_hbonds)
}