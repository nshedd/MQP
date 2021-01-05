marker_elements = read.table(path.expand("~/pfc_marker_peaks.txt"), sep = '\t', header=TRUE, row.names=1)
head(marker_elements)

top_10_elements = marker_elements[1,]

for (i in 1:12) {
  temp_elements = marker_elements[marker_elements$group == i,]
  temp_elements_sorted = temp_elements[order(temp_elements$Log2FC, decreasing=TRUE),]
  j = nrow(temp_elements_sorted)
  if (j >= 3) {
    top_temp_elements = temp_elements_sorted[1:3,]
    top_10_elements = rbind(top_10_elements, top_temp_elements)
  }
  else if (j>0) {
    top_temp_elements = temp_elements_sorted[1:j,]
    top_10_elements = rbind(top_10_elements, top_temp_elements)
  } 
}

top_10_elements = top_10_elements[-1,]

write.table(top_10_elements, file = path.expand("~/pfc_marker_peaks_top3.txt"), sep = '\t')
