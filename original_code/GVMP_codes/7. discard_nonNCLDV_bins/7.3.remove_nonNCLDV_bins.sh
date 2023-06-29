# input file : 'NonNCLDV_bins_list.txt'
## Remove non-NCLDV bins fna files
#!/bin/bash

while IFS= read -r bin_id; do
    # Add the '.fna' suffix to the bin_id
    filename="${bin_id}.fna" #You can replace with '.faa' if there are any faa files of the non-NCLDV bins you want to remove as well
    
    # Check if the file exists and remove it
    if [ -f "$filename" ]; then
        rm "$filename" #you can add -f flag if you want to remove them for sure
        echo "Removed file: $filename"
    fi
done < NonNCLDV_bins_list.txt
