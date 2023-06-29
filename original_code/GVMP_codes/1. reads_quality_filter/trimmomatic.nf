#input file:*_1|2.fastq.gz
#output file: *_1|2.paired.fastq.gz/ *_1|2.unpaired.fastq.gz

params.query = "/some/data/*_1.fastq.gz"
params.outputpath = "/some/data/"

reads = file(params.query)

process trimmomatic_batch {

   [ directives ]

   input:
    path(reads)

   output:
    < process outputs >

   script:
   """
   module load trimmoatic
   find ${query}| parallel --dry-run "java -jar /usr/appli/freeware/trimmomatic/0.38/trimmomatic-0.38.jar PE -threads 6 {} {}_2 /outputpath{/.}.paired.fastq.gz /outputpath/{/.}.unpaired.fastq.gz /outputpath/{/.}_2.paired.fastq.gz /outputpath/{/.}_2.unpaired.fastq.gz ILLUMINACLIP:/usr/appli/freeware/trimmomatic/0.38/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36” > trimmomatic.qsubarray
   """
}