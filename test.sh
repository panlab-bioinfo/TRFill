declare -i end_dis
chrs=( chr1, chr1, chr1, chr2)
gap_starts=(100 500 1000 2000)
gap_ends=(200 700 1500 3000)
end_dis=500
i=0
for ((j=$((i+1)); j<${#chrs[@]}; j++))
do
    if [ ${chrs[$j]} == ${chrs[$i]} ]; then
        gap_starts[$j]=$((gap_starts[$j] + end_dis))
        gap_ends[$j]=$((gap_ends[$j] + end_dis))
    else
        break
    fi
done
echo ${gap_starts[@]}
echo ${gap_ends[@]}