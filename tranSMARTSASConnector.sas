/******************************************************/
/* SAS Connector for Transmart through ODBC or Oracle */
/*         ConvergeHEALTH, David Newton 2014          */
/******************************************************/

%macro tm_listStudies(out=studylist, search= , conceptsize=4, gex=false, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let sql = (SELECT DISTINCT SOURCESYSTEM_CD as STUDYCODE, SUBSTR(C_FULLNAME,1,REGEXP_INSTR(C_FULLNAME,'\\',1,&conceptsize.)) as CONCEPT_PATH FROM I2B2 WHERE SOURCESYSTEM_CD LIKE (&quote.&wildcard.&search.&wildcard.&quote.) ORDER BY SOURCESYSTEM_CD);
    
    %tm_query(out=&out, sql=&sql, debug=&debug);

	%IF &gex=true %THEN %DO;
		%let dscount = 0;
		%tm_dataSetSize(var=dscount, dataset=&out);
		%IF &dscount = 0 %THEN
			%DO;
				%RETURN;
			%END;
		%ELSE
			%DO;
				%let css = ;
				%tm_datasetToDelimitedString(var=css, column=STUDYCODE, dataset=&out, delimiter='&quote.,&quote.', distinct=true);
				%let css = &quote.&css.&quote.;
			%END;

		%let sql = (SELECT COUNT (DISTINCT PROBESET_ID) AS PROBE_COUNT, TRIAL_NAME as STUDYCODE FROM deapp.de_subject_microarray_data WHERE TRIAL_NAME IN (&css.) GROUP BY TRIAL_NAME ORDER BY TRIAL_NAME);
        
        %tm_query(out=work.probeCounts, sql=&sql, debug=&debug);

		data &out;
			merge &out work.probecounts;
			by STUDYCODE;
		run;
	%END;
%mend tm_listStudies;

%macro tm_listSearchTerms(out=searchtermlist, search= , category=, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let search = %qupcase(&search);
	%let category = %qupcase(&category);
	%let sql = SELECT DISTINCT skt.keyword_Term, skt.rank, sk.data_category FROM SEARCH_KEYWORD_TERM skt join SEARCHAPP.search_keyword sk on skt.search_keyword_id = sk.search_keyword_id WHERE skt.keyword_Term LIKE(&quote.&search.&wildcard.&quote.);
	/* Note no wildcard at start of search term */
	%if %sysevalf(%superq(category)=,boolean) %then
		%do;
		%end;
	%else
		%do;
			%let sql = &sql AND sk.data_Category = &quote.&category.&quote.;
		%end;
	%if &debug=true %then %put &sql;
    
    %tm_query(out=&out, sql=(&sql), debug=&debug);

%mend tm_listSearchTerms;

%macro tm_getProbeGeneSNPMapping(out=pgsnpmapping, probeIds=null, geneIds=null, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let sql = (SELECT sgm.snp_name SNP, bm.BIO_MARKER_NAME GENE FROM de_snp_gene_map sgm INNER JOIN bio_marker bm ON bm.PRIMARY_EXTERNAL_ID = to_char(sgm.ENTREZ_GENE_ID) WHERE ;

	%if %sysfunc(exist(&probeIds)) & %sysfunc(exist(&geneIds)) %then
		%do;
			%put tm_getProbeGeneSNPMapping: Please provide ONLY either a list of probe IDs or gene IDs.;
			%return;
		%end;
	
	%if %sysfunc(exist(&probeIds)) %then
		%do;
			%let pidcss = ;
			%tm_datasetToDelimitedString(var=pidcss, dataset=&probeIds, delimiter='&quote.,&quote.');
			%let pidcss = &quote.&pidcss.&quote.;
			%put &pidcss;
			%let sql = &sql sgm.snp_name IN (&pidcss));
		%end;
	%else %if %sysfunc(exist(&geneIds)) %then
		%do;
			%let gidcss = ;
			%tm_datasetToDelimitedString(var=gidcss, dataset=&geneIds, delimiter='&quote.,&quote.');
			%let gidcss = %qupcase(&quote.&gidcss.&quote.);
			%put &gidcss;
			%let sql = &sql bm.bio_marker_name IN (&gidcss));
		%end;
	%else
		%do;
			%put tm_getProbeGeneSNPMapping: Please provide either a list of probe IDs or gene IDs.;
			%return;
		%end;
        
    %tm_query(out=&out, sql=&sql, debug=&debug);
    
%mend tm_getProbeGeneSNPMapping;

%macro tm_getProbeGeneMapping(out=pgmapping, probeIds=null, geneIds=null, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let sql = (SELECT DISTINCT PROBE_ID, GENE_SYMBOL FROM de_mrna_annotation WHERE ;

	%if %sysfunc(exist(&probeIds)) & %sysfunc(exist(&geneIds)) %then
		%do;
			%put tm_getProbeGeneMapping: Please provide ONLY either a list of probe IDs or gene IDs.;
			%return;
		%end;
	
	%if %sysfunc(exist(&probeIds)) %then
		%do;
			%let pidcss = ;
			%tm_datasetToDelimitedString(var=pidcss, dataset=&probeIds, delimiter='&quote.,&quote.');
			%let pidcss = &quote.&pidcss.&quote.;
			%put &pidcss;
			%let sql = &sql probe_id IN (&pidcss));
		%end;
	%else %if %sysfunc(exist(&geneIds)) %then
		%do;
			%let gidcss = ;
			%tm_datasetToDelimitedString(var=gidcss, dataset=&geneIds, delimiter='&quote.,&quote.');
			%let gidcss = %qupcase(&quote.&gidcss.&quote.);
			%put &gidcss;
			%let sql = &sql gene_symbol IN (&gidcss));
		%end;
	%else
		%do;
			%put tm_getProbeGeneMapping: Please provide either a list of probe IDs or gene IDs.;
			%return;
		%end;

    %tm_query(out=&out, sql=&sql, debug=&debug);
    
%mend tm_getProbeGeneMapping;

%macro tm_getDistinctConcepts(out=distinctConcepts, studyList=null, pathMatchList=null, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let sql = (SELECT CD.CONCEPT_CD, CD.CONCEPT_PATH, CD.SOURCESYSTEM_CD as STUDYCODE, COUNT(1) as COUNT FROM CONCEPT_DIMENSION CD INNER JOIN OBSERVATION_FACT OBSF ON OBSF.CONCEPT_CD = CD.CONCEPT_CD WHERE ;

	%if %sysfunc(exist(&studyList)) %then
		%do;
			%let studycss = ;
			%tm_datasetToDelimitedString(var=studycss, dataset=&studyList, delimiter='&quote.,&quote.');
			%let studycss = %qupcase(&quote.&studycss.&quote.);
			%if &debug=true %then %put tm_getDistinctConcepts: Study string: &studycss;
			%let sql = &sql CD.SOURCESYSTEM_CD IN (&studycss) AND ;
		%end;

	%if ~%sysfunc(exist(&pathMatchList)) %then
		%do;
			%put tm_getDistinctConcepts: Please provide a pathMatchList.;
			%return;
		%end;
		
	%let pathcss = ;
	%tm_datasetToDelimitedString(var=pathcss, dataset=&pathMatchList, delimiter='&wildcard.&quote. OR upper(CD.CONCEPT_PATH) LIKE &quote.&wildcard.');
	%let pathcss = %qupcase(&quote.&wildcard.&pathcss.&wildcard.&quote.);
	%if &debug=true %then %put tm_getDistinctConcepts: Path match string: &pathcss;

	%let sql = &sql (;
	%let sql = &sql upper(CD.CONCEPT_PATH) LIKE &pathcss;
	%let sql = &sql );
		
	%let sql = &sql GROUP BY CD.CONCEPT_PATH, CD.CONCEPT_CD, CD.SOURCESYSTEM_CD);


	%put &sql;

    %tm_query(out=&out, sql=&sql, debug=&debug);
%mend tm_getDistinctConcepts;

%macro tm_getClinicalData(out=clinicalData, conceptCodeList=, pivot=true, prePivotTrim=true, trimLength=4, debug=false) / store;
	%put In getClinicalData;
	%let wildcard = %;
	%let quote = %str(%');
	%let sql = (SELECT DISTINCT PT.PATIENT_NUM, REPLACE(PD.SOURCESYSTEM_CD,PT.TRIAL || ':') SUBJECT_ID, PT.TRIAL TRIAL_NAME, CD.CONCEPT_PATH, CD.CONCEPT_CD, CD.NAME_CHAR, case OBSF.VALTYPE_CD WHEN 'T' THEN TVAL_CHAR WHEN 'N' THEN CAST(NVAL_NUM AS varchar2(30)) END VALUE FROM CONCEPT_DIMENSION CD INNER JOIN OBSERVATION_FACT OBSF ON OBSF.CONCEPT_CD = CD.CONCEPT_CD INNER JOIN PATIENT_TRIAL PT  ON PT.PATIENT_NUM = OBSF.PATIENT_NUM INNER JOIN PATIENT_DIMENSION PD ON PD.PATIENT_NUM = PT.PATIENT_NUM ;
    
	%let ccCount = 0;
	%tm_dataSetSize(dataset=&conceptCodeList, var=ccCount);
	%if &ccCount=0 %then 
		%do;
			%if &debug=true %then %put tm_getClinicalData: No concept codes provided, will return a blank data set;
			%let sql = &sql WHERE CD.CONCEPT_CD IN (''));
		%end;
	%else
		%do;
			%let conceptcss = ;
			%tm_datasetToDelimitedString(var=conceptcss, dataset=&conceptCodeList, delimiter='&quote.,&quote.');
			%let conceptcss = &quote.&conceptcss.&quote.;
			%if &debug=true %then %put &conceptcss;
			%let sql = &sql WHERE CD.CONCEPT_CD IN (&conceptcss));
		%end;

    %tm_query(out=&out, sql=&sql, debug=&debug);

	/* We have the data - pivot if necessary */
	%if &pivot=true %then
		%do;
			/* Take off the last \folder\ of the concept path */
			%if &prePivotTrim=true %then
				%do;
					data &out(drop=regex);
						set &out;
						CONCEPT_PATH=prxchange('s/\\[^\\]+\\?$//', -1, trim(CONCEPT_PATH));
						/* If trim length is above 0, take folders off the start */
						%if &trimLength>0 %then
							%do;
								regex = cat('s/^\\(.*?\\){', &trimLength, '}//');
								CONCEPT_PATH=prxchange(regex, -1, trim(CONCEPT_PATH));	
							%end;
						/* If less than 0, take folders off the end */
						%else %if &trimLength<0 %then
							%let trimLength = 0-&trimLength;
							%do;
								regex = cat('s/(\\([^\\])*){', &trimLength, '}$//');
								CONCEPT_PATH=substr(trim(CONCEPT_PATH), prxmatch(regex, trim(CONCEPT_PATH))+1);	
							%end;
					run;
				%end;

			proc sort data=&out out=&out;
				by PATIENT_NUM;
			run;

			proc transpose data=&out out=&out;
				var VALUE;
				by PATIENT_NUM SUBJECT_ID TRIAL_NAME;
				id CONCEPT_PATH;
			run;
		%end;
%mend tm_getClinicalData;

%macro tm_getPatientMapping(out=patientMapping, studyList=, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let sql = (SELECT PT.TRIAL TRIAL_NAME, REPLACE(PD.SOURCESYSTEM_CD,PT.TRIAL || ':') SUBJECT_ID, PD.PATIENT_NUM PATIENT_ID FROM PATIENT_TRIAL PT INNER JOIN PATIENT_DIMENSION PD ON PD.PATIENT_NUM = PT.PATIENT_NUM ;
	%let studycss = ;
	%tm_datasetToDelimitedString(var=studycss, dataset=&studyList, delimiter='&quote.,&quote.');
	%let studycss = %qupcase(&quote.&studycss.&quote.);
	%let sql = &sql WHERE PT.TRIAL IN (&studycss));
	
    %tm_query(out=&out, sql=&sql, debug=&debug);
%mend tm_getPatientMapping;

%macro tm_listHDDAttributes(sampleTypes=sampleTypes, tissueTypes=tissueTypes, timepoints=timepoints, studyList=, debug=) / store;
	%let wildcard = %;
	%let quote = %str(%');

	%let studycss = ;
	%tm_datasetToDelimitedString(var=studycss, dataset=&studyList, delimiter='&quote.,&quote.');
	%let studycss = %qupcase(&quote.&studycss.&quote.);

	%tm_listHDDAttributesInternal(out=&sampleTypes,studycss=&studycss,sqlColumn=SAMPLE_TYPE,debug=&debug);
	%tm_listHDDAttributesInternal(out=&tissueTypes,studycss=&studycss,sqlColumn=TISSUE_TYPE,debug=&debug);
	%tm_listHDDAttributesInternal(out=&timepoints,studycss=&studycss,sqlColumn=TIMEPOINT,debug=&debug);

%mend tm_listHDDAttributes;

%macro tm_getSNPData(out=snpData, studyList=null, geneList=null, pathway=null, signature=null, patientList=null, sampleTypes=null, tissueTypes=null, timepoints=null, platforms=null, probeList=null, showGenes=false, pivot=true, cnPivotAggregate=null, gtPivotAggregate=null, pivotPatientId=false, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');

	%if ~%sysfunc(exist(&studyList)) & ~%sysfunc(exist(&patientList)) %then
		%do;
			%put tm_getSNPData: Please provide either a list of studies or a list of patients.;
			%return;
		%end;
        
	%if ~%sysfunc(exist(&pathway)) & ~%sysfunc(exist(&signature)) & ~%sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getSNPData: Please provide a pathway, signature or gene list.;
			%return;
		%end;

	%if %sysfunc(exist(&pathway)) & %sysfunc(exist(&signature)) %then
		%do;
			%put tm_getSNPData: You cannot filter by both a pathway and a gene signature.;
			%return;
		%end;

	%if %sysfunc(exist(&pathway)) & %sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getSNPData: You cannot filter by both a pathway and a gene list.;
			%return;
		%end;

	%if %sysfunc(exist(&signature)) & %sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getSNPData: You cannot filter by both a signature and a gene list.;
			%return;
		%end;

	%let baseSqlSelect = SNP_GENO.SNP_NAME AS SNP, DSM.PATIENT_ID, DSM.SUBJECT_ID, DSM.sample_type,DSM.timepoint,DSM.tissue_type,SNP_GENO.SNP_CALLS AS GENOTYPE,SNP_COPY.COPY_NUMBER AS COPYNUMBER,PD.sourcesystem_cd,DSM.GPL_ID;
	%let baseSqlTable = FROM DE_SUBJECT_SAMPLE_MAPPING DSM INNER JOIN patient_dimension PD ON DSM.patient_id = PD.patient_num LEFT JOIN DE_SNP_CALLS_BY_GSM SNP_GENO ON DSM.OMIC_PATIENT_ID = SNP_GENO.PATIENT_NUM AND DSM.SAMPLE_CD = SNP_GENO.GSM_NUM LEFT JOIN DE_SNP_COPY_NUMBER SNP_COPY ON DSM.OMIC_PATIENT_ID = SNP_COPY.PATIENT_NUM AND SNP_GENO.snp_name = SNP_COPY.snp_name LEFT JOIN DE_SNP_PROBE PR ON PROBE_NAME=SNP_GENO.snp_name ;
	%let baseSqlWhere = ;
	
	/* Filter by patient or study list, whichever isn't blank */
	%if %sysfunc(exist(&studyList)) %then
		%do;
			%let studycss = ;
			%tm_datasetToDelimitedString(var=studycss, dataset=&studyList, delimiter='&quote.,&quote.');
			%let studycss = %qupcase(&quote.&studycss.&quote.);
			%let baseSqlWhere = WHERE DSM.trial_name IN (&studycss);
		%end;
	%else
		%do;
			%if ~%sysfunc(exist(&patientList)) %then
				%do;
					%put tm_getSNPData: You cannot filter by both a signature and a gene list.;
					%return;
				%end;
			%else
				%do;
					%let patientcss = ;
					%tm_datasetToDelimitedString(var=patientcss, dataset=&patientList, delimiter='&quote.,&quote.');
					%let patientcss = %qupcase(&quote.&patientcss.&quote.);
					%let baseSqlWhere = WHERE DSM.patient_id IN (&patientcss);
				%end;
		%end;
	
	/* Append more conditions if a pathway was included */
	%if %sysfunc(exist(&pathway)) or %sysfunc(exist(&geneList)) or %sysfunc(exist(&probeList)) %then
		%do;
			%let baseSqlSelect = SELECT &baseSqlSelect;
			%if &showGenes=true %then %let baseSqlSelect = &baseSqlSelect, bm.BIO_MARKER_NAME AS GENE ;
			%let newTableJoins = ;
			%if %sysfunc(exist(&pathway)) or %sysfunc(exist(&geneList)) %then
				%do;
					%let newTableJoins = INNER JOIN DE_SNP_GENE_MAP D2 ON D2.SNP_NAME = PR.SNP_NAME INNER JOIN bio_marker bm ON bm.PRIMARY_EXTERNAL_ID = to_char(D2.ENTREZ_GENE_ID) INNER JOIN bio_marker_correl_mv sbm ON sbm.asso_bio_marker_id = bm.bio_marker_id INNER JOIN search_keyword sk ON sk.bio_data_id = sbm.bio_marker_id INNER JOIN SEARCH_KEYWORD_TERM skt ON sk.SEARCH_KEYWORD_ID = skt.SEARCH_KEYWORD_ID ;
				%end;
			
			/* Create the string of genes (pathway or list) or probes */
			%let searchWords = ;
			%if %sysfunc(exist(&pathway)) %then
				%do;
					%tm_datasetToDelimitedString(var=searchWords, dataset=&pathway, delimiter='&quote.,&quote.');
					%let searchWords = %qupcase(&quote.&searchWords.&quote.);
					%let baseSqlWhere = &baseSqlWhere AND skt.KEYWORD_TERM IN (&searchWords);
				%end;
			%else %if %sysfunc(exist(&geneList)) %then
				%do;
					%tm_datasetToDelimitedString(var=searchWords, dataset=&geneList, delimiter='&quote.,&quote.');
					%let searchWords = %qupcase(&quote.&searchWords.&quote.);
					%let baseSqlWhere = &baseSqlWhere AND skt.KEYWORD_TERM IN (&searchWords);
				%end;
			%else %if %sysfunc(exist(&probeList)) %then
				%do;
					%tm_datasetToDelimitedString(var=searchWords, dataset=&probeList, delimiter='&quote.,&quote.');
					%let searchWords = &quote.&searchWords.&quote.; /* No uppercasing here! */
					%let baseSqlWhere = &baseSqlWhere AND SNP_GENO.SNP_NAME IN (&searchWords);
				%end;
			
			%let baseSqlTable = &baseSqlTable &newTableJoins;
		%end;
	%else %if %sysfunc(exist(&signature)) %then
		%do;
			%let baseSqlSelect = SELECT &baseSqlSelect;
			%if &showGenes=true %then %let baseSqlSelect = &baseSqlSelect, bm.BIO_MARKER_NAME AS GENE ;
			%let newTableJoins = INNER JOIN DE_SNP_GENE_MAP D2 ON D2.SNP_NAME = PR.SNP_NAME INNER JOIN bio_marker bm ON bm.PRIMARY_EXTERNAL_ID = to_char(D2.ENTREZ_GENE_ID) INNER JOIN search_bio_mkr_correl_fast_mv sbm ON sbm.asso_bio_marker_id = bm.bio_marker_id INNER JOIN search_keyword sk ON sk.bio_data_id = sbm.domain_object_idINNER JOIN SEARCH_KEYWORD_TERM skt ON sk.SEARCH_KEYWORD_ID = skt.SEARCH_KEYWORD_ID ;
			%let baseSqlTable = &baseSqlTable &newTableJoins;
			%let searchWords = ;
			%tm_datasetToDelimitedString(var=searchWords, dataset=&signature, delimiter='&quote.,&quote.');
			%let searchWords = %qupcase(&quote.&searchWords.&quote.);
			%let baseSqlWhere = &baseSqlWhere AND skt.KEYWORD_TERM IN (&searchWords);
		%end;
	%else
		%do;
			%let baseSqlSelect = SELECT &baseSqlSelect;
			%if &showGenes=true %then %let baseSqlSelect = &baseSqlSelect, bm.BIO_MARKER_NAME AS GENE ;
			%let newTableJoins = INNER JOIN DE_SNP_GENE_MAP D2 ON D2.SNP_NAME = PR.SNP_NAME INNER JOIN bio_marker bm ON bm.PRIMARY_EXTERNAL_ID = to_char(D2.ENTREZ_GENE_ID) ;
			%let baseSqlTable = &baseSqlTable &newTableJoins;
		%end;

	%if %sysfunc(exist(&sampleTypes)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&sampleTypes, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND DSM.sample_type IN (&css);
		%end;

	%if %sysfunc(exist(&tissueTypes)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&tissueTypes, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND DSM.tissue_type IN (&css);
		%end;

	%if %sysfunc(exist(&timepoints)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&timepoints, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND DSM.timepoint IN (&css);
		%end;

	%if %sysfunc(exist(&platforms)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&platforms, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND DSM.GPL_ID IN (&css);
		%end;
	
	/* At last, put together the final query! */

	%let sql = (&baseSqlSelect &baseSqlTable &baseSqlWhere);

    %tm_query(out=&out, sql=&sql, debug=&debug);

	/* Now pivot the data set, if required */
	%if &pivot=false %then %return;

	%let pivotPatient = SUBJECT_ID;
	%if &pivotPatientId=true %then %let pivotPatient = PATIENT_ID;

	%let columnsCN = &pivotPatient, COPYNUMBER, SNP, GPL_ID;
	%let columnsGT = &pivotPatient, GENOTYPE, SNP, GPL_ID;
	
	%let pivotFields = GPL_ID SNP;
	%let sortField = GPL_ID;

	%if showGenes=true %then
		%do;
			%let columnsCN = &columnsCN., GENE;
			%let columnsGT = &columnsGT., GENE;
			%let pivotFields = GENE GPL_ID SNP;
			%let sortField = GENE;
		%end;

	/* Call the provided aggregate/pivot functions, if both are provided */
	%if ~(&cnPivotAggregate=null and &gtPivotAggregate=null) %then
		%do;
			%if &cnPivotAggregate=null or &gtPivotAggregate=null %then
				%do;
					%put tm_getSNPData: If specifying a function to aggregate and pivot, please specify for both genotype (gtPivotAggregate) and copy number (cnPivotAggregate).;
					%return;
				%end;
				%if &debug=true %then %put tm_getSNPData: Using passed-in macros to aggregate and pivot the data.;
                
				%&cnPivotAggregate.(out=&out.CN);
				%&gtPivotAggregate.(out=&out.GT);
				%return;
		%end;

	%if &debug=true %then %put tm_getSNPData: Using plain pivot with no aggregate.;

	/* If no aggregate functions provided, just pivot. */

    proc sql;
        create table &out.CN as
        select &columnsCN. from &out;
    quit;

    proc sql;
        create table &out.GT as
        select &columnsGT. from &out;
    quit;
                
            
	proc sort data=&out.CN out=&out.CN;
		by SNP;
	run;

	proc transpose data=&out.CN out=&out.CN;
		var COPYNUMBER;
		by &pivotFields.;
		id &pivotPatient.;
	run;

	proc sort data=&out.GT out=&out.GT;
		by SNP;
	run;

	proc transpose data=&out.GT out=&out.GT;
		var GENOTYPE;
		by &pivotFields.;
		id &pivotPatient.;
	run;

%mend tm_getSNPData;

%macro tm_getGEXData(out=gexData, studyList=null, geneList=null, pathway=null, signature=null, patientList=null, sampleTypes=null, tissueTypes=null, timepoints=null, platforms=null, probeList=null, removeOverlappingPlatforms=null, showGenes=false, pivot=true, pivotAggregate=null, pivotPatientId=false, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');

	%tm_getPartitionList(studyList=&studyList, debug=&debug);

	%if ~%sysfunc(exist(&studyList)) & ~%sysfunc(exist(&patientList)) %then
		%do;
			%put tm_getSNPData: Please provide either a list of studies or a list of patients.;
			%return;
		%end;
        
	%if ~%sysfunc(exist(&pathway)) & ~%sysfunc(exist(&signature)) & ~%sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getSNPData: Please provide a pathway, signature or gene list.;
			%return;
		%end;

	%if %sysfunc(exist(&pathway)) & %sysfunc(exist(&signature)) %then
		%do;
			%put tm_getSNPData: You cannot filter by both a pathway and a gene signature.;
			%return;
		%end;

	%if %sysfunc(exist(&pathway)) & %sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getSNPData: You cannot filter by both a pathway and a gene list.;
			%return;
		%end;

	%if %sysfunc(exist(&signature)) & %sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getSNPData: You cannot filter by both a signature and a gene list.;
			%return;
		%end;

	%let baseSqlSelect = a.PATIENT_ID, a.RAW_INTENSITY, a.ZSCORE, a.LOG_INTENSITY, a.assay_id, ssm.subject_id, ssm.sample_type, ssm.timepoint, ssm.tissue_type, ssm.trial_name, ssm.GPL_ID, b.probe_id,  b.probeset_id;
	%let baseSqlTable = FROM deapp.de_subject_microarray_data a INNER JOIN de_subject_sample_mapping ssm ON ssm.assay_id = A.assay_id;
	%let baseSqlWhere = ;
	
	/* Filter by patient or study list, whichever isn't blank */
	%if %sysfunc(exist(&studyList)) %then
		%do;
			%let studycss = ;
			%tm_datasetToDelimitedString(var=studycss, dataset=&studyList, delimiter='&quote.,&quote.');
			%let studycss = %qupcase(&quote.&studycss.&quote.);
			%let studyIndexcss = ;
			%tm_datasetToDelimitedString(var=studyIndexcss, dataset=partitionList, delimiter='&quote.,&quote.');
			%let studyIndexcss = &quote.&studyIndexcss.&quote.;
			%let baseSqlWhere = WHERE SSM.trial_name IN (&studycss) AND a.trial_source in (&studyIndexcss);
		%end;
	%else
		%do;
			%let patientcss = ;
			%tm_datasetToDelimitedString(var=patientcss, dataset=&patientList, delimiter='&quote.,&quote.');
			%let patientcss = %qupcase(&quote.&patientcss.&quote.);
			%let baseSqlWhere = WHERE SSM.patient_id IN (&patientcss);
		%end;
	
	/* Append more conditions if a pathway was included */
	%if %sysfunc(exist(&pathway)) or %sysfunc(exist(&geneList)) or %sysfunc(exist(&probeList)) %then
		%do;
			%let baseSqlSelect = SELECT DISTINCT &baseSqlSelect;
			%if &showGenes=true %then %let baseSqlSelect = &baseSqlSelect, b.GENE_SYMBOL, b.GENE_ID ;
			%let newTableJoins = ;
			%if %sysfunc(exist(&pathway)) or %sysfunc(exist(&geneList)) %then
				%do;
					%let newTableJoins = INNER JOIN deapp.de_mrna_annotation b ON a.probeset_id = b.probeset_id INNER JOIN bio_marker bm ON bm.PRIMARY_EXTERNAL_ID = to_char(b.GENE_ID) INNER JOIN bio_marker_correl_mv sbm ON sbm.asso_bio_marker_id = bm.bio_marker_id INNER JOIN search_keyword sk ON sk.bio_data_id = sbm.bio_marker_id INNER JOIN SEARCH_KEYWORD_TERM skt ON sk.SEARCH_KEYWORD_ID = skt.SEARCH_KEYWORD_ID ;
					%let baseSqlTable = &baseSqlTable &newTableJoins;
				%end;
			
			/* Create the string of genes (pathway or list) or probes */
			%let searchWords = ;
			%if %sysfunc(exist(&pathway)) %then
				%do;
					%tm_datasetToDelimitedString(var=searchWords, dataset=&pathway, delimiter='&quote.,&quote.');
					%let searchWords = %qupcase(&quote.&searchWords.&quote.);
					%let baseSqlWhere = &baseSqlWhere AND skt.KEYWORD_TERM IN (&searchWords);
				%end;
			%else %if %sysfunc(exist(&geneList)) %then
				%do;
					%tm_datasetToDelimitedString(var=searchWords, dataset=&geneList, delimiter='&quote.,&quote.');
					%let searchWords = %qupcase(&quote.&searchWords.&quote.);
					%let baseSqlWhere = &baseSqlWhere AND skt.KEYWORD_TERM IN (&searchWords);
				%end;
			%else %if %sysfunc(exist(&probeList)) %then
				%do;
					%tm_datasetToDelimitedString(var=searchWords, dataset=&probeList, delimiter='&quote.,&quote.');
					%let searchWords = &quote.&searchWords.&quote.; /* No uppercasing here! */
					%let baseSqlWhere = &baseSqlWhere AND b.PROBE_ID IN (&searchWords);
				%end;
			
			
		%end;
	%else %if %sysfunc(exist(&signature)) %then
		%do;
			%let baseSqlSelect = SELECT DISTINCT &baseSqlSelect;
			%if &showGenes=true %then %let baseSqlSelect = &baseSqlSelect, b.GENE_SYMBOL, b.GENE_ID ;
			%let newTableJoins = INNER JOIN (SELECT CASE WHEN bfg.PROBESET_ID IS NOT NULL THEN bfg.PROBESET_ID WHEN bbm.PROBESET_ID IS NOT NULL THEN bbm.PROBESET_ID END PROBESET_ID, ;
			%let newTableJoins = &newTableJoins CASE WHEN bfg.PROBESET_ID IS NOT NULL THEN bfg.PROBE_ID WHEN bbm.PROBESET_ID IS NOT NULL THEN bbm.PROBE_ID END PROBE_ID, ;
			%let newTableJoins = &newTableJoins CASE WHEN bfg.PROBESET_ID IS NOT NULL THEN bfg.GENE_SYMBOL WHEN bbm.PROBESET_ID IS NOT NULL THEN bbm.GENE_SYMBOL END GENE_SYMBOL ;
			%let newTableJoins = &newTableJoins FROM SEARCH_KEYWORD_TERM skt INNER JOIN search_keyword sk ON sk.SEARCH_KEYWORD_ID = skt.SEARCH_KEYWORD_ID INNER JOIN SEARCH_GENE_SIGNATURE SGS ON SGS.SEARCH_GENE_SIGNATURE_ID = sk.bio_data_id INNER JOIN SEARCH_GENE_SIGNATURE_ITEM SGSI ON SGS.SEARCH_GENE_SIGNATURE_ID = SGSI.SEARCH_GENE_SIGNATURE_ID ;
			%let newTableJoins = &newTableJoins LEFT JOIN bio_assay_feature_group fg ON fg.bio_assay_feature_group_id = SGSI.bio_assay_feature_group_id LEFT JOIN deapp.de_mrna_annotation bfg ON bfg.PROBE_ID = fg.FEATURE_GROUP_NAME LEFT JOIN bio_marker bm ON bm.bio_marker_id = SGSI.bio_marker_id LEFT JOIN deapp.de_mrna_annotation bbm ON to_char (bbm.GENE_ID) = bm.PRIMARY_EXTERNAL_ID ;

			%let searchWords = ;
			%tm_datasetToDelimitedString(var=searchWords, dataset=&signature, delimiter='&quote.,&quote.');
			%let searchWords = %qupcase(&quote.&searchWords.&quote.);

			%let newTableJoins = &newTableJoins WHERE skt.KEYWORD_TERM IN (&signature)) b ON a.PROBESET_ID = b.PROBESET_ID ;

			%let baseSqlTable = &baseSqlTable &newTableJoins;
		%end;
	%else
		%do;
			%let baseSqlSelect = SELECT DISTINCT &baseSqlSelect;
			%let newTableJoins = INNER JOIN deapp.de_mrna_annotation b ON a.probeset_id = b.probeset_id ;
			%let baseSqlTable = &baseSqlTable &newTableJoins;
		%end;

	%if %sysfunc(exist(&sampleTypes)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&sampleTypes, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND SSM.sample_type IN (&css);
		%end;

	%if %sysfunc(exist(&tissueTypes)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&tissueTypes, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND SSM.tissue_type IN (&css);
		%end;

	%if %sysfunc(exist(&timepoints)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&timepoints, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND SSM.timepoint IN (&css);
		%end;

	%if %sysfunc(exist(&platforms)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&platforms, delimiter='&quote.,&quote.');
			%let css = &quote.&css.&quote.;
			%let baseSqlWhere = &baseSqlWhere AND SSM.GPL_ID IN (&css);
		%end;
	
	/* At last, put together the final query! */

	%let sql = (&baseSqlSelect &baseSqlTable &baseSqlWhere);

    %tm_query(out=&out, sql=&sql, debug=&debug);

	/* Now pivot the data set, if required */
	%if &pivot=true %then 
		%do;

			%let pivotPatient = SUBJECT_ID;
			%if &pivotPatientId=true %then %let pivotPatient = PATIENT_ID;

			%let columnsGEX = &pivotPatient, LOG_INTENSITY, PROBE_ID, GPL_ID;
			%let pivotFields = GPL_ID PROBE_ID;
			%let sortFields = GPL_ID, PROBE_ID;

			%if &showGenes=true %then
				%do;
					%let columnsGEX = &columnsGEX., GENE_SYMBOL;
					%let pivotFields = GENE_SYMBOL GPL_ID PROBE_ID;
					%let sortFields = GENE_SYMBOL, GPL_ID, PROBE_ID;
				%end;

			/* Call the provided aggregate/pivot functions, if both are provided */
			%if &pivotAggregate=null %then
				%do;
					%if &debug=true %then %put tm_getGEXData: Using plain pivot with no aggregate.;

					proc sql;
						create table &out as
						select &columnsGEX from &out ORDER BY &sortFields;
					quit;

					proc transpose data=&out out=&out;
						var LOG_INTENSITY;
						by &pivotFields.;
						id &pivotPatient.;
					run;
				%end;
			%else
				%do;
					%if &debug=true %then %put tm_getGEXData: Using passed-in macro to aggregate and pivot the data.;
					%&pivotAggregate.(out=&out);
					%return;
				%end;
		%end;

	%if %sysfunc(exist(&removeOverlappingPlatforms)) %then
		%do;
			/* Get a list of all probes that occur more than once in the data set */
			proc sql;
				create table probeFrequency as
				select PROBE_ID as Name, COUNT(*) AS FREQ FROM &out GROUP BY PROBE_ID;
			quit;

			proc sql;
				create table probeFrequency as
				select * from probeFrequency WHERE FREQ > 1;
			quit;
			
			%let probecss = ;
			%tm_datasetToDelimitedString(var=probecss, dataset=probeFrequency, delimiter='","');
			%let probecss = "&probecss.";

			%let platformcss = ;
			%tm_datasetToDelimitedString(var=platformcss, dataset=&removeOverlappingPlatforms, delimiter='","');
			%let platformcss = "&platformcss.";

			%let sql = select * from &out WHERE PROBE_ID NOT IN (&probecss) OR GPL_ID NOT IN (&platformcss);
			%if &debug=true %then %put sql;

			proc sql;
				create table &out as
				(&sql);
			quit;
		%end;

%mend tm_getGEXData;

%macro tm_getGeneGoMembership(out=geneGoMembership, geneGoName=null, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	
	%let baseSqlSelect = SELECT DISTINCT BMA.BIO_MARKER_NAME AS PATHWAY_NAME, BMB.BIO_MARKER_NAME AS GENE;
	%let baseSqlTable = FROM BIO_MARKER BMA, BIO_DATA_CORRELATION BDC, BIO_MARKER BMB;
	%let baseSqlWhere = WHERE BMA.BIO_MARKER_ID = BDC.BIO_DATA_ID and BDC.ASSO_BIO_DATA_ID = BMB.BIO_MARKER_ID and BMA.PRIMARY_SOURCE_CODE = 'GO';
	%let baseSqlSort = order by BMA.BIO_MARKER_NAME, BMB.BIO_MARKER_NAME;

	%if %sysfunc(exist(&geneGoName)) %then
		%do;
			%let css = ;
			%tm_datasetToDelimitedString(var=css, dataset=&geneGoName, delimiter='&quote.,&quote.');
			%let css = %qupcase(&quote.&css.&quote.);
			%let baseSqlWhere = &baseSqlWhere AND UPPER(BMA.bio_marker_name) IN (&css);
		%end;

	%let sql = (&baseSqlSelect &baseSqlTable &baseSqlWhere &baseSqlSort);

    %tm_query(out=&out, sql=&sql, debug=&debug);
%mend tm_getGeneGoMembership;

%macro tm_getClinicalMutationData(out=clinicalMutationData, studyList=null, geneList=null, trimLength=4, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let gene =;

	%if ~%sysfunc(exist(&studyList)) %then
		%do;
			%put tm_getClinicalMutationData: Please provide a studyList.;
			%return;
		%end;

	%if ~%sysfunc(exist(&geneList)) %then
		%do;
			%put tm_getClinicalMutationData: Please provide a geneList.;
			%return;
		%end;

	%if &debug=true %then %put tm_getClinicalMutationData: Getting patient mapping for studies;
	%tm_getPatientMapping(out=work.mutationPatientMapping, studyList=&studyList, debug=&debug);
	
	proc sort data=work.mutationPatientMapping out=work.mutationPatientMapping;
		by PATIENT_ID;
	run;

	%let patientCount = 0;
	%tm_dataSetSize(dataset=work.mutationPatientMapping, var=patientCount);
	%if &patientCount=0 %then 
		%do;
			%put tm_getClinicalMutationData: No patients found for the given study list: &studyList;
			%return;
		%end;

	proc sql noprint;
		select Name into :genecss separated by '||'
		from &geneList;
	quit;

	%let geneLength=;
	%tm_dataSetSize(var=geneLength, dataset=&geneList);

	%let i=1;
	%do %while(&i <= &geneLength); /* Indexes from 1 */
		%let currentGene = %scan(&genecss.,&i);
		%if &debug=true %then %put tm_getClinicalMutationData: Getting clinical mutation data for &currentGene;
		%tm_getClinicalMutDataPerGene(out=&out, studyList=&studyList, trimLength=&trimLength, gene=&currentGene, debug=&debug);
		%if &debug=true %then %put tm_getClinicalMutationData: Merging clinical mutation data for &currentGene;
		/*Merge in clinical mutation data*/
		data work.mutationPatientMapping;
			merge work.mutationPatientMapping work.mutationDataUnPivotSort;
			by PATIENT_ID;
		run;

		%if &debug=true %then %put tm_getClinicalMutationData: Getting allele frequency data for &currentGene;
		%tm_getClinicalAllFreqPerGene(out=&out, studyList=&studyList, trimLength=&trimLength, gene=&currentGene, debug=&debug);
		%if &debug=true %then %put tm_getClinicalMutationData: Merging allele frequency data for &currentGene;
		data work.mutationPatientMapping;
			merge work.mutationPatientMapping work.alleleFrequencyUnPivotSort;
			by PATIENT_ID;
		run;

		%let i = &i + 1;
	%end;

	proc sql noprint;
		create table &out as
		select * from work.mutationPatientMapping;
	quit;

%mend tm_getClinicalMutationData;

/* UTILITY MACROS */

%macro tm_datasetToDelimitedString(var=css, column=, dataset=, delimiter=',', distinct=false) / store;
		
		%if %sysevalf(%superq(dataset)=,boolean) %then
			%do;
				%put No data set provided to convert to string;
			%return;
			%end;

		/* If no column provided, use the first. */
		%if %sysevalf(%superq(column)=,boolean) %then
			%do;
				proc sql noprint;
				data _null_;
				  set &dataset (obs=1);
				  array c{*} _CHARACTER_;
				  call symput('column',vname(c{1}));
				run;
			%end;
		%if %sysevalf(%superq(column)=,boolean) %then
			%do;
				proc sql noprint;
				data _null_;
				  set &dataset (obs=1);
				  array n{*} _NUMERIC_;
				  call symput('column',vname(n{1}));
				run;
			%end;

        %if &distinct=true %then
            %do;
                proc sql noprint;
                    select distinct &column
                    into :&var
                    separated by &delimiter
                    from &dataset;
                quit;
            %end;
        %else
            %do;
                proc sql noprint;
                    select &column
                    into :&var
                    separated by &delimiter
                    from &dataset;
                quit;
            %end;
%mend tm_datasetToDelimitedString;

%macro tm_dataSetSize(var=, dataset=) / store;
		proc sql noprint;
			select count(*)
			into :&var
			from &dataset;
		quit;
%mend tm_dataSetSize;

%macro tm_listHDDAttributesInternal(out=,studycss=,sqlColumn=,debug=) / store;
	%let sql = (SELECT DISTINCT &sqlColumn FROM DE_SUBJECT_SAMPLE_MAPPING WHERE TRIAL_NAME IN (&studycss));
    %tm_query(out=&out, sql=&sql, debug=&debug);
%mend tm_listHDDAttributesInternal;

%macro tm_getClinicalMutDataPerGene(out=geneMutationData, studyList=, trimLength=4, gene=, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let parameterList = Mutations&wildcard.&gene.&wildcard.Mutation Type;

	data work.mutationDataPerGenePath;
		length Name $32;
		Name = "&parameterList";
		output;
	run;

	%tm_getDistinctConcepts(out=mutationDataPerGeneConcepts, studyList=&studyList, pathMatchList=work.mutationDataPerGenePath, debug=&debug);

	data mutationDataPerGeneConcepts;
		set mutationDataPerGeneConcepts;
		rename CONCEPT_CD=Name;
		output;
	run;

	%tm_getClinicalData(out=mutationDataUnPivot, pivot=false, conceptCodeList=mutationDataPerGeneConcepts, debug=&debug);

	data mutationDataUnPivot(drop=regex);
		set mutationDataUnPivot;
		&gene=trim(CONCEPT_PATH);
		/* If trim length is above 0, take folders off the start */
		%if &trimLength>0 %then
			%do;
				regex = cat('s/^\\(.*?\\){', &trimLength, '}//');
				&gene=prxchange(regex, -1, &gene);	
			%end;
		/* If less than 0, take folders off the end */
		%else %if &trimLength<0 %then
			%let trimLength = 0-&trimLength;
			%do;
				regex = cat('s/(\\([^\\])*){', &trimLength, '}$//');
				&gene=substr(&gene, prxmatch(regex, &gene)+1);	
			%end;

		PATIENT_ID=PATIENT_NUM;
	run;

	proc sql noprint;
		CREATE TABLE mutationDataUnPivotSort AS
		select PATIENT_ID,&gene
		from mutationDataUnPivot
		order by PATIENT_ID;
	quit;

%mend tm_getClinicalMutDataPerGene;

%macro tm_getClinicalAllFreqPerGene(out=alleleFrequencyData, studyList=, trimLength=4, gene=, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');
	%let parameterList = Mutations&wildcard.&gene.&wildcard.Mutant Allele Frequency;

	data work.alleleFrequencyPerGenePath;
		length Name $32;
		Name = "&parameterList";
		output;
	run;

	%tm_getDistinctConcepts(out=alleleFrequencyPerGeneConcepts, studyList=&studyList, pathMatchList=work.alleleFrequencyPerGenePath, debug=&debug);

	data alleleFrequencyPerGeneConcepts;
		set alleleFrequencyPerGeneConcepts;
		rename CONCEPT_CD=Name;
		output;
	run;

	%tm_getClinicalData(out=alleleFrequencyUnPivot, pivot=false, conceptCodeList=alleleFrequencyPerGeneConcepts, debug=&debug);

	data alleleFrequencyUnPivot(drop=regex);
		set alleleFrequencyUnPivot;
		&gene._Allele_Freq=trim(CONCEPT_PATH);
		/* If trim length is above 0, take folders off the start */
		%if &trimLength>0 %then
			%do;
				regex = cat('s/^\\(.*?\\){', &trimLength, '}//');
				&gene._Allele_Freq=prxchange(regex, -1, &gene._Allele_Freq);	
			%end;
		/* If less than 0, take folders off the end */
		%else %if &trimLength<0 %then
			%let trimLength = 0-&trimLength;
			%do;
				regex = cat('s/(\\([^\\])*){', &trimLength, '}$//');
				&gene._Allele_Freq=substr(&gene._Allele_Freq, prxmatch(regex, &gene._Allele_Freq)+1);	
			%end;

		PATIENT_ID=PATIENT_NUM;
	run;

	proc sql noprint;
		CREATE TABLE alleleFrequencyUnPivotSort AS
		select PATIENT_ID,&gene._Allele_Freq
		from alleleFrequencyUnPivot
		order by PATIENT_ID;
	quit;

%mend tm_getClinicalAllFreqPerGene;

%macro tm_getPartitionList(studyList=, debug=false) / store;
	%let wildcard = %;
	%let quote = %str(%');

	%let partsql = SELECT DISTINCT omic_source_study||&quote.:&quote.||source_cd as Name from de_subject_sample_mapping WHERE trial_name in ;
	%let studycss = ;
	%tm_datasetToDelimitedString(var=studycss, dataset=&studyList, delimiter='&quote.,&quote.');
	%let studycss = %qupcase(&quote.&studycss.&quote.);
	%let partsql = (&partsql(&studycss));

    %tm_query(out=partitionList, sql=&partsql, debug=&debug);

%mend tm_getPartitionList;

/*
This macro is an example for connecting to an ODBC data source called TSMART. Uncomment it and remove
the Oracle function below if you want to use ODBC.

%macro tm_query(out=out, sql=, debug=false) / store;

    %if &debug=true %then %put TM running SQL: &sql;
	proc sql;
		CONNECT TO ODBC AS TSMART(database=&tmdsn user=&tmuser password=&tmpassword);
		CREATE TABLE &out AS
		SELECT * FROM CONNECTION TO TSMART
		&sql;
	quit;
%mend;
*/

%macro tm_query(out=out, sql=, debug=false) / store;
    %if &debug=true %then %put TM running SQL: &sql;
    proc sql;
        CONNECT TO ORACLE as ORACLE(user=&tmuser password=&tmpassword path=&tmdsn);
        CREATE TABLE &out AS
        SELECT * FROM CONNECTION TO ORACLE
        &sql;
    quit;
%mend;
