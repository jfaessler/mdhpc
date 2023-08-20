for size in {1..16}
do
  echo "Processing cluster with $size layers"
  formatted=$(printf "%02d" "$size")
  ./07 generated_clusters/cluster_"$formatted"l.xyz data_"$formatted"l.csv > out_"$formatted"l.out
done