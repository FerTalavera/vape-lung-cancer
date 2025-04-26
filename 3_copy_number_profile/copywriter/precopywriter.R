library('CopyhelpeR')
library('CopywriteR')

preCopywriteR(output.folder = file.path("/mnt/atgc-d1/drobles/ftalavera/vape_lung_cancer/3_copy_number_profile/copywriter/precopywriter"),
                bin.size = 20000,
                ref.genome = "mm10",
                prefix = "")
