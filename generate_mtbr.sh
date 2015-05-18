#!/bin/bash

usage()
{
	echo "usage: `basename $0` mtbr.conf in.sam"
}

checkerror()
{
	progress=$1
	
	if [ -z "`grep -i error ${log_file}`" ] && [ -z "`grep -i "No such file or directory" ${log_file}`" ]
	then
		echo "${progress}=OK" >> ${progress_file} 		
	else
		echo "ERROR : current progress \"${progress}\""
		exit -1
	fi
}

sub_track()
{		
	strand=$1
	track_file=$2
	if [ "${strand}" == "p" ];then
		full_strand="positive"
	else
		full_strand="negative"
	fi
	
		
	echo  "" >> ${track_file}
    echo -e "\t\ttrack ${group_info}.${full_strand}" >> ${track_file}
	echo -e "\t\tbigDataUrl ${bigwig_path}/${track_group}_${strand}.bw" >> ${track_file}
	echo -e "\t\tshortLabel ${group_info}.${strand}" >> ${track_file}
	echo -e "\t\tlongLabel ${group_info}.${full_strand}" >> ${track_file}
	echo -e "\t\tparent ${group_info}-separated" >> ${track_file}
	echo -e "\t\ttype bigWig" >> ${track_file}
	echo -e "\t\tcolor 0,0,255" >> ${track_file}
	echo -e "\t\taltColor 0,0,255" >> ${track_file}
}

gen_trackhub()
{
		
	track_data_path=${trackhub_path}"/"${genome_type}
	mkdir -p ${track_data_path}
	track_bigwig_path=${track_data_path}"/"${bigwig_path}
	rm -rf ${track_bigwig_path}
	mv ${bigwig_path} ${track_data_path}"/" 

	# hub.txt 

	hub_file=${trackhub_path}"/hub.txt"
	echo "hub hub_name(change this to your projectname)" > ${hub_file}
	echo "shortLabel hub_short_label" >> ${hub_file}
	echo "longLabel hub_long_label" >> ${hub_file}
	echo "genomesFile genomes.txt" >> ${hub_file}
	echo "email ${USER}_lilab@icsc.dlmedu.edu.cn" >> ${hub_file}
	echo "descriptionUrl hub.html" >> ${hub_file}

	# genomes.txt

	genomes_file=${trackhub_path}"/genomes.txt"
	echo "genome ${genome_type}" > ${genomes_file}
	echo "trackDb ${genome_type}/trackDb.txt" >> ${genomes_file}

	# hub.html

	hub_html=${trackhub_path}"/hub.html"
	echo "<!DOCTYPE html>" > ${hub_html}
	echo "<html>" >> ${hub_html}
	echo "<head>" >> ${hub_html}
	echo "Track Hub Description" >> ${hub_html}
	echo "</head>" >> ${hub_html}
	echo "<body>" >> ${hub_html}
	echo "<p> Your track hub description here. </p>" >> ${hub_html}
	echo "</body>" >> ${hub_html}
	echo "</html>" >> ${hub_html}

	# trackDb.txt
	trackdb_txt=${track_data_path}"/trackDb.txt"
	track_groups=`find ${track_bigwig_path} -maxdepth 1 -type f -regextype posix-extended -regex ".*\_[np]\.bw" | xargs -I {} basename {} | cut -d "_" -f1 | sort | uniq`
	
	for track_group in ${track_groups};do
		group_info=${sam_file_base}"."${track_group}
		data_type=${track_group##*.}
		
		echo ""	>> ${trackdb_txt}
		echo "####################################" >> ${trackdb_txt}
		echo ""	>> ${trackdb_txt}
		echo "track ${group_info}-separated"  >> ${trackdb_txt}
		echo "type bigWig" >> ${trackdb_txt}
		echo "container multiWig" >> ${trackdb_txt}
		echo "shortLabel ${group_info}.separated" >> ${trackdb_txt}
		echo "longLabel ${group_info}.separated" >> ${trackdb_txt}
		echo "visibility hide" >> ${trackdb_txt}
		echo "aggregate transparentOverlay" >> ${trackdb_txt}
		echo "showSubtrackColorOnUi on" >> ${trackdb_txt}
		echo "yLineOnOff on" >> ${trackdb_txt}
		echo "priority 1.0" >> ${trackdb_txt}
	
		if [ "${data_type}" == "score" ]
		then
			echo "viewLimits -1:1"  >> ${trackdb_txt}
			echo "autoScale off" >> ${trackdb_txt}
		else 
			echo "autoScale on" >> ${trackdb_txt}
		fi
		
		echo "configurable on" >> ${trackdb_txt}
		echo "windowingFunction mean" >> ${trackdb_txt}
		
		sub_track "p" ${trackdb_txt}
		sub_track "n" ${trackdb_txt}		
	done
	
	track_data_list=`find ${track_bigwig_path} -maxdepth 1 -type f -regextype posix-extended -regex "[^\_]*\.bw" | xargs -I {} basename {}`
	
	for track_data in ${track_data_list};do
		group_info=${sam_file_base}"."${track_data}
		data_type=`echo ${track_data} | rev | cut -d"." -f-2 | rev | cut -d "." -f1`
		
		echo ""	>> ${trackdb_txt}
		echo "####################################" >> ${trackdb_txt}
		echo ""	>> ${trackdb_txt}
		echo "track ${group_info}-mixed"  >> ${trackdb_txt}
		echo "bigDataUrl ${bigwig_path}/${track_data}" >> ${trackdb_txt}
		echo "type bigWig" >> ${trackdb_txt}
		echo "shortLabel ${group_info}.mixed" >> ${trackdb_txt}
		echo "longLabel ${group_info}.mixed" >> ${trackdb_txt}
		echo "visibility hide" >> ${trackdb_txt}
		echo "priority 1.0" >> ${trackdb_txt}
		echo "yLineOnOff on" >> ${trackdb_txt}
		echo "color 255,0,0" >> ${trackdb_txt}
		
		if [ "${data_type}" == "score" ]
		then
			echo "viewLimits 0:1"  >> ${trackdb_txt}
			echo "autoScale off" >> ${trackdb_txt}
		else 
			echo "autoScale on" >> ${trackdb_txt}
		fi
		
		echo "configurable on" >> ${trackdb_txt}
		echo "windowingFunction mean" >> ${trackdb_txt}
	done

}

if [ $# != 2 ] 
then
	usage
	exit -1
fi

script_path=`dirname $0`
conf_file=$1
sam_file=$2

callmeth="perl ${script_path}/utils/callMajorMethyl.pl"
samtoolscmd="samtools" 
cgmtbr_r="${script_path}/Rscript/get_cg_mtbr.Rscript"
wigtobigwig="${script_path}/utils/wigToBigWig"

# chmod tools 
chmod +x ${cgmtbr_r}
chmod +x ${wigtobigwig}

# soure the config file

if [ ! -f ${conf_file} ]
then
	echo "[error]: \"${conf_file}\" no such file!"
	exit -1
fi
source ${conf_file}

# check samtools version

samtools_version=`${samtoolscmd} --version 2>&1 | grep samtools | cut -d- -f1 | cut -d" " -f2`

if [ -z "${samtools_version}" ] || [ `expr ${samtools_version} \< 1.1` != 0 ]
then
	echo "samtools version must great than 1.2"
	exit -1
fi

# check fa file

if [ ! -f ${ref_fa} ]
then 
	echo " You must set a valide ref file in ${conf_file} "
	echo " error ref file is ${ref_fa} "
	exit -1
fi

# check get_cg_mtbr.Rscript options

cgmtbr_r_option=""
if [ -n "${genome_library}" ] && [ -n "${genome_name}"  ]
then
	cgmtbr_r_option=" -l ${genome_library} -n ${genome_name} "
else
	if [ -n "${genome_type}" ]
	then
		cgmtbr_r_option=" -t ${genome_type} "
	else
		echo "genome_type or genome_library & genome_name must be set in ${conf_file}"
		echo "error genome_type is ${genome_type}"
		echo "error genome_library is ${genome_library}"
		echo "error genome_name is ${genome_name}"
		exit -1
	fi
fi

cgmtbr_r_option=${cgmtbr_r_option}" -a ${position_offset} -d ${mtbr_header} -w ${wig_type} "

# check input sam file

if [ ! -f ${sam_file} ]
then
	echo "[error]: sam file \"${sam_file}\" does not exist!"
	exit -1
fi

# generate c2t & g2a files & paths

sam_file_base=`basename ${sam_file}`
sam_file_base=${sam_file_base%.*}

temp_ct_sam_file="tmp."${sam_file_base}".ct.sam"
temp_ga_sam_file="tmp."${sam_file_base}".ga.sam"
temp_ct_bam_file="tmp."${sam_file_base}".ct.bam"
temp_ga_bam_file="tmp."${sam_file_base}".ga.bam"
temp_ct_sorted_file="tmp."${sam_file_base}".ct.sorted.bam"
temp_ga_sorted_file="tmp."${sam_file_base}".ga.sorted.bam"
temp_ct_pileup_file="tmp."${sam_file_base}".ct.pileup"
temp_ga_pileup_file="tmp."${sam_file_base}".ga.pileup"


log_file=${sam_file_base}".log"
mtbr_path=${sam_file_base}".mtbr.all"

# check and load progress file

progress_file="."${sam_file_base}".progress"

if [ -f ${progress_file} ]
then
	source ${progress_file}
else
	touch ${progress_file}
fi

# do the conversions

echo "[*] Processing `basename ${sam_file}` `date`" | tee ${log_file}

echo "[*] Generating C2T & G2A sam files..." | tee -a ${log_file}

if [ -z ${progress_split_bsmark} ]
then
	${samtoolscmd} view -hS ${sam_file} | awk '/^@/{print > "'${temp_ct_sam_file}'"; print > "'${temp_ga_sam_file}'";} /XB:Z:F.*CT/{ print >> "'${temp_ct_sam_file}'";}/XB:Z:F.*GA/{print >> "'${temp_ga_sam_file}'" }' >> ${log_file} 2>&1
	checkerror "progress_split_bsmark"
fi

echo "[*] Converting & Sorting C2T & G2A bam files..." | tee -a ${log_file}

if [ -z ${progress_sort_sam} ]
then
	${samtoolscmd} view -b ${temp_ct_sam_file} -o ${temp_ct_bam_file} >> ${log_file} 2>&1
	${samtoolscmd} sort -m 5G ${temp_ct_bam_file} -o ${temp_ct_sorted_file} -T ${temp_ct_sorted_file} >> ${log_file} 2>&1
	${samtoolscmd} view -b ${temp_ga_sam_file} -o ${temp_ga_bam_file} >> ${log_file} 2>&1
	${samtoolscmd} sort -m 5G ${temp_ga_bam_file} -o ${temp_ga_sorted_file} -T ${temp_ga_sorted_file} >> ${log_file} 2>&1
	echo "progress_sort_sam=OK" >> ${progress_file}
fi

echo "[*] Generating MTBRs..." | tee -a ${log_file}

if [ -z ${progress_generate_mtbr} ]
then

	rm -rf ${mtbr_path}
	mkdir ${mtbr_path}

	# pipup with chromfa
	temp_meth_ct_file="tmp."${sam_file_base}".ct.meth" 
	temp_meth_ga_file="tmp."${sam_file_base}".ga.meth" 

	echo "    pileup reads with chromfas..." | tee -a ${log_file}
	${samtoolscmd} mpileup -B -f ${ref_fa} ${temp_ct_sorted_file} -o ${temp_ct_pileup_file} >> ${log_file} 2>&1
	${callmeth} -s 1 -i ${temp_ct_pileup_file} -o ${temp_meth_ct_file} >> ${log_file} 2>&1

	${samtoolscmd} mpileup -B -f ${ref_fa} ${temp_ga_sorted_file} -o ${temp_ga_pileup_file} >> ${log_file} 2>&1
	${callmeth} -s 2 -i ${temp_ga_pileup_file} -o ${temp_meth_ga_file} >> ${log_file} 2>&1
	
	checkerror "progress_generate_mtbr"
fi

# generate mtbr
echo "[*] Spliting mtbr file..." | tee -a ${log_file}

if [ -z ${progress_split_mtbr} ]	
then
	cd ${mtbr_path}
	awk 'BEGIN{ chr = "" }
		{
			if (chr != $1){ 
				chr = $1; 
				print $0 >> chr;	 	
			}else{
				print $0 >> chr;
			}
		}' "../"${temp_meth_ct_file} "../"${temp_meth_ga_file}
	cd ".."
	checkerror "progress_split_mtbr" 
fi

echo "[*] Generating cg mtbr file and wig file" | tee -a ${log_file} 


# generate Wig file
if [ -z ${progress_generate_wig} ]	
then
	${cgmtbr_r} ${mtbr_path} ${cgmtbr_r_option}  >> ${log_file} 2>&1 
	checkerror "progress_generate_wig"
fi


# Build trackhub
echo "[*] Buliding trackhub files..." | tee -a ${log_file}

wig_path=${sam_file_base}".mtbr.wig"
bigwig_path=${sam_file_base}".mtbr.bigwig"
trackhub_path=${sam_file_base}".trackhub"

mkdir ${bigwig_path}

if [ -z ${progress_build_trackhub} ]
then
	if [ -z "${genome_type}" ]
	then

		while [ "${genome_type}" != "mm9"  ] && [ "${genome_type}" != "mm10" ] && [ "${genome_type}" != "hg18" ] && [ "${genome_type}" != "hg19" ] 
		do
			read -p "Require genome type [ mm9 | mm10 | hg18 | hg19 ]: " genome_type
		done
	fi
	
	chrom_sizes="${script_path}/chrom.sizes/${genome_type}.chrom.sizes"	
	wig_files=`find ${wig_path} -maxdepth 1 -type f -name "*.wig"`
	
	for wig_file in ${wig_files};do
		echo "    wig to bigwig : ${wig_file}" | tee -a ${log_file}
		wig_base=`basename ${wig_file}`
		wig_base=${wig_base%.*}
		bigwig_file=${bigwig_path}"/"${wig_base}".bw"
		${wigtobigwig} ${wig_file} ${chrom_sizes} ${bigwig_file}  >> ${log_file} 2>&1  
	done
		
	gen_trackhub 
	checkerror "progress_generate_wig"
fi




# delete temp files
rm -f tmp*
rm ${progress_file}

echo "[*] Complete `date`" | tee -a ${log_file}
