package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.utilities.workflows.OicrWorkflow;
import java.util.Map;
import java.util.logging.Logger;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;

/**
 * <p>
 * For more information on developing workflows, see the documentation at
 * <a href="http://seqware.github.io/docs/6-pipeline/java-workflows/">SeqWare
 * Java Workflows</a>.</p>
 *
 * Quick reference for the order of methods called: 1. setupDirectory 2.
 * setupFiles 3. setupWorkflow 4. setupEnvironment 5. buildWorkflow
 *
 * See the SeqWare API for
 * <a href="http://seqware.github.io/javadoc/stable/apidocs/net/sourceforge/seqware/pipeline/workflowV2/AbstractWorkflowDataModel.html#setupDirectory%28%29">AbstractWorkflowDataModel</a>
 * for more information.
 */
public class Varscan2Workflow extends OicrWorkflow {

    //dir
    private String dataDir, tmpDir;
    private String outDir;

    // Input Data
    private String tumorBam;
    private String normalBam;
    private String outputFilenamePrefix;

    //varscan intermediate file names
    private String snpFile;
    private String cnvFile;
    private String indelFile;
    private String mpileupFile;
    private String varscanCopycallFile;
    private String copyNumberFile;
    private String somaticPileupFile;

    //Tools
    private String samtools;
    private String java;
    private String varscan;

    //Memory allocation
    private Integer varscanMem;
    private String javaMem = "-Xmx8g";

    //path to bin
    private String bin;
//    private String pypy;
    private String rScript;
    private String rLib;

    //ref Data
    private String refFasta;
    private String intervalFile;
//    private String sequenzaGCData;

    private boolean manualOutput;
    private static final Logger logger = Logger.getLogger(Varscan2Workflow.class.getName());
    private String queue;
    private Map<String, SqwFile> tempFiles;

    // meta-types
    private final static String TXT_METATYPE = "text/plain";
    private final static String TAR_GZ_METATYPE = "application/tar-gzip";

    private void init() {
        try {
            //dir
            dataDir = "data";
            tmpDir = getProperty("tmp_dir");

            // input samples 
            tumorBam = getProperty("input_files_tumor");
            normalBam = getProperty("input_files_normal");

            //Ext id
            outputFilenamePrefix = getProperty("external_name");

            //samtools
            samtools = getProperty("samtools");
            java = getProperty("java");
            varscan = getProperty("varscan").toString();

            // ref fasta
            refFasta = getProperty("ref_fasta");
            intervalFile = getProperty("interval_bed");

            manualOutput = Boolean.parseBoolean(getProperty("manual_output"));
            queue = getOptionalProperty("queue", "");

            varscanMem = Integer.parseInt(getProperty("varscan_mem"));

        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void setupDirectory() {
        init();
        this.addDirectory(dataDir);
        this.addDirectory(tmpDir);
        if (!dataDir.endsWith("/")) {
            dataDir += "/";
        }
        if (!tmpDir.endsWith("/")) {
            tmpDir += "/";
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {
        SqwFile file0 = this.createFile("tumor");
        file0.setSourcePath(tumorBam);
        file0.setType("application/bam");
        file0.setIsInput(true);
        SqwFile file1 = this.createFile("normal");
        file1.setSourcePath(normalBam);
        file1.setType("application/bam");
        file1.setIsInput(true);
        return this.getFiles();
    }

    @Override
    public void buildWorkflow() {

        /**
         * Steps for sequenza: 1. Check if "bam" file exists; true 2. Check if
         * "bai" file exists; true: go to step 4 3. Check if normal Pb_R sample
         * exists; true: go to step 4; else abort 3. If false: samtools index
         * "bam" file 4. Run job sequenza-utils 5. If outputFile ends with
         * "bin50.gz"; go to step 6; else go to step 4 6. Run job sequenzaR 7.
         * Iterate through the files/folders in outDir: 8. If fileName1 ==
         * "pandc.txt" and fileName2 ends with "Total_CN.seg"; create a folder
         * called "copynumber" 9. If fileType == "folder"; create a folder
         * called "model-fit"; move folders to "model-fit" 10. If fileType ==
         * "file" && fileName != outputFile; move file to "model-fit" 11. Delete
         * outputFile (rm outputFile) 12. zip "model-fit" 13. outputFile =
         * fileName2 14. OutputDir contains the following: fileName1,
         * outputFile, model-fit.zip
         */
        // workflow : read inputs tumor and normal bam; run sequenza-utils; write the output to temp directory; 
        // run sequenzaR; handle output; provision files (3) -- model-fit.zip; text/plain; text/plain
        Job parentJob = null;
        this.outDir = this.outputFilenamePrefix + "_output";
        this.snpFile = this.tmpDir + this.outputFilenamePrefix + ".varscanSomatic.snp";
        this.cnvFile = this.tmpDir + this.outputFilenamePrefix + ".VarScan.CopyNumber.copynumber";
        this.indelFile = this.tmpDir + this.outputFilenamePrefix + ".varscanSomatic.indel";
        this.mpileupFile = this.tmpDir + this.outputFilenamePrefix + ".mpileup";
        this.varscanCopycallFile = this.tmpDir + this.outputFilenamePrefix + ".VarScan.CopyCaller";
        this.copyNumberFile = this.tmpDir + this.outputFilenamePrefix + ".VarScan.CopyNumber";
        this.somaticPileupFile = this.tmpDir + this.outputFilenamePrefix + ".varscanSomatic";

        Job mpileup = runMpileup();
        parentJob = mpileup;

        Job somaticMpileup = getSomaticPileup();
        somaticMpileup.addParent(parentJob);
        parentJob = somaticMpileup;

        Job varscanIndels = varscanIndels();
        varscanIndels.addParent(parentJob);

        Job varscanSNP = varscanSNP();
        varscanSNP.addParent(parentJob);

        Job varscanCNA = varscanCNA();
        varscanCNA.addParent(parentJob);
        parentJob = varscanCNA;

        Job varscanCNACaller = varscanCNACaller();
        varscanCNACaller.addParent(parentJob);
        parentJob = varscanCNACaller;

        // Provision all files
        SqwFile snpOutFile = createOutputFile(this.snpFile, TXT_METATYPE, this.manualOutput);
        snpOutFile.getAnnotations().put("SNPS from VarScan2 ", "Varscan_SNPS ");
        varscanSNP.addFile(snpOutFile);
        SqwFile indelsFile = createOutputFile(this.indelFile, TXT_METATYPE, this.manualOutput);
        indelsFile.getAnnotations().put("Indels from VarScan2 ", "Varscan_Indels ");
        varscanIndels.addFile(indelsFile);
        SqwFile copyNumverVarFile = createOutputFile(this.cnvFile, TXT_METATYPE, this.manualOutput);
        copyNumverVarFile.getAnnotations().put("CNA from VarScan2 ", "Varscan_CNA ");
        varscanCNA.addFile(copyNumverVarFile);
        SqwFile copyNumFile = createOutputFile(this.copyNumberFile, TXT_METATYPE, this.manualOutput);
        copyNumFile.getAnnotations().put("CNA calls from VarScan2 ", "Varscan_CNA_calls ");
        varscanCNA.addFile(copyNumFile);
        SqwFile copyCallFile = createOutputFile(this.copyNumberFile, TXT_METATYPE, this.manualOutput);
        copyCallFile.getAnnotations().put("CNV calls from VarScan2 ", "Varscan_copy_number_calls ");
        varscanCNACaller.addFile(copyCallFile);
    }

    
    private Job runMpileup() {
        Job mpileup = getWorkflow().createBashJob("mpileup");
        Command cmd = mpileup.getCommand();
        cmd.addArgument(this.samtools);
        cmd.addArgument("mpileup -q 1");
        cmd.addArgument("-f " + this.refFasta);
        cmd.addArgument("-l " + this.intervalFile);
        cmd.addArgument("-B -d 1000000");
        cmd.addArgument(getFiles().get("normal").getProvisionedPath());
        cmd.addArgument(getFiles().get("tumor").getProvisionedPath());
        cmd.addArgument("> " + this.mpileupFile);
        mpileup.setMaxMemory(Integer.toString(varscanMem * 1024));
        mpileup.setQueue(getOptionalProperty("queue", ""));
        return mpileup;
    }

    private Job getSomaticPileup() {
        Job somaticPileup = getWorkflow().createBashJob("somatic_pileup");
        Command cmd = somaticPileup.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar");
        cmd.addArgument(this.varscan);
        cmd.addArgument("somatic");
        cmd.addArgument(this.mpileupFile);
        cmd.addArgument(this.somaticPileupFile);
        cmd.addArgument("--mpileup 1");
        somaticPileup.setMaxMemory(Integer.toString(varscanMem * 1024));
        somaticPileup.setQueue(getOptionalProperty("queue", ""));
        return somaticPileup;
    }

    private Job varscanIndels() {
        Job vsIndels = getWorkflow().createBashJob("varscan_indels");
        Command cmd = vsIndels.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar " + this.varscan);
        cmd.addArgument("processSomatic");
        cmd.addArgument(this.indelFile);
        vsIndels.setMaxMemory(Integer.toString(varscanMem * 1024));
        vsIndels.setQueue(getOptionalProperty("queue", ""));
        return vsIndels;
    }

    private Job varscanSNP() {
        Job vsSNP = getWorkflow().createBashJob("varscan_snps");
        Command cmd = vsSNP.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar " + this.varscan);
        cmd.addArgument("processSomatic");
        cmd.addArgument(this.snpFile);
        vsSNP.setMaxMemory(Integer.toString(varscanMem * 1024));
        vsSNP.setQueue(getOptionalProperty("queue", ""));
        return vsSNP;
    }

    private Job varscanCNA() {
        Job vsCNA = getWorkflow().createBashJob("varscan_cna");
        Command cmd = vsCNA.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar " + this.varscan);
        cmd.addArgument("copynumber");
        cmd.addArgument(this.mpileupFile);
        cmd.addArgument(this.copyNumberFile);
        cmd.addArgument("--mpileup 1");
        vsCNA.setMaxMemory(Integer.toString(varscanMem * 1024));
        vsCNA.setQueue(getOptionalProperty("queue", ""));
        return vsCNA;
    }

    private Job varscanCNACaller() {
        Job vsCNAcall = getWorkflow().createBashJob("varscan_cna_call");
        Command cmd = vsCNAcall.getCommand();
        cmd.addArgument(this.java);
        cmd.addArgument(this.javaMem);
        cmd.addArgument("-jar " + this.varscan);
        cmd.addArgument("copyCaller");
        cmd.addArgument(this.cnvFile);
        cmd.addArgument("--output-file " + this.varscanCopycallFile);
        cmd.addArgument("--mpileup 1");
        vsCNAcall.setMaxMemory(Integer.toString(varscanMem * 1024));
        vsCNAcall.setQueue(getOptionalProperty("queue", ""));
        return vsCNAcall;
    }


}
