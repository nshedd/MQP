marker_elements = read.table(path.expand("~/temporal_marker_peaks.txt"), sep = '\t', header=TRUE, row.names=1)
head(marker_elements)

top_10_elements = marker_elements[1,]

for (i in 1:10) {
  temp_elements = marker_elements[marker_elements$group == i,]
  print(temp_elements)
  temp_elements_sorted = order(temp_elements$Log2FC)
  j = nrow(temp_elements_sorted)
  print(j)
  if (j >= 10) {
    top_temp_elements = temp_elements_sorted[1:10,]
  }
  else {
    top_temp_elements = temp_elements_sorted[1:j,]
  }
  top_10_elements = rbind(top_10_elements, top_temp_elements)
}

top_10_elements = top_10_elements[!1,]

write.table(markerList, file = path.expand("~/temporal_marker_peaks_top10.txt"), sep = '\t')
