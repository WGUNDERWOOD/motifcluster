# Compress pdf files

for filename in R/doc/*.pdf
do
    echo Compressing $filename ...
    optpdf $filename
    optpdf $filename
    optpdf $filename
done

for filename in R/vignettes/*.pdf
do
    echo Compressing $filename ...
    optpdf $filename
    optpdf $filename
    optpdf $filename
done

for filename in python/tutorial/*.pdf
do
    echo Compressing $filename ...
    optpdf $filename
    optpdf $filename
    optpdf $filename
done
