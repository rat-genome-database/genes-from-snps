package edu.mcw.rgd.genesFromSnps;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.datamodel.Eva;
import edu.mcw.rgd.datamodel.RgdId;
import edu.mcw.rgd.datamodel.SpeciesType;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.datamodel.variants.VariantSampleDetail;
import edu.mcw.rgd.process.Utils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.beans.factory.support.DefaultListableBeanFactory;
import org.springframework.beans.factory.xml.XmlBeanDefinitionReader;
import org.springframework.core.io.FileSystemResource;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class Main {
    private String version;
    protected Logger logger = LogManager.getLogger("status");
    private DAO dao = new DAO();
    private SimpleDateFormat sdt = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
    private Integer sampleId;

    public static void main(String[] args) throws Exception {
        DefaultListableBeanFactory bf = new DefaultListableBeanFactory();
        new XmlBeanDefinitionReader(bf).loadBeanDefinitions(new FileSystemResource("properties/AppConfigure.xml"));
        try {
            Main main = (Main) bf.getBean("main");
            for (int i = 0; i < args.length; i++) {
                switch (args[i]) {
                    default:
                        main.run();
                        return;
                    case "-getRsIds":
                        main.run2();
                        return;
                    case "-close":
                        main.closestGene();
                        return;

                }
            }
        }
        catch (Exception e) {
            Utils.printStackTrace(e, LogManager.getLogger("status"));
            throw e;
        }
    }

    void run() throws Exception {
        File folder = new File("data/human_snps/");
        List<File> files = listFilesInFolder(folder);
        geneCacheMap = new HashMap<>();
        logger.info(version);
        long pipeStart = System.currentTimeMillis();
        logger.info("Pipeline started at "+sdt.format(new Date(pipeStart))+"\n");
//        System.out.println(files.size());
        for (File file : files){
            BufferedReader br = openFile(file.getAbsolutePath());
            logger.info("\tGetting genes for: "+file.getName());
            BufferedWriter bw = new BufferedWriter(new FileWriter("genes_from_"+file.getName()));
            String lineData;
            while ((lineData = br.readLine()) != null) {
                // get gene in snp region
                if (lineData.startsWith("Locus")) {
                    bw.write(lineData+"\tGene Closest\tGenes Upstream\tGenes Downstream\n");
                    continue;
                }
                bw.write(lineData);
                String[] cols = lineData.split("\t");
                String chr = cols[1];
                int pos = Integer.parseInt(cols[2]);
                int boundLimit = 0;
                List<Integer> closeGeneIds = getClosestGene(boundLimit, pos, chr);
                while (closeGeneIds.isEmpty()) {
                    boundLimit = boundLimit + 1000;
                    closeGeneIds = getClosestGene(boundLimit, pos, chr);
                }

                String closeGeneList = listOfGenesToPrint(closeGeneIds);
                bw.write("\t"+closeGeneList);
//                if (!geneIds.isEmpty()){
//                    logger.debug("\t\tGetting Gene");
//                    String geneList = listOfGenesToPrint(geneIds);
//                    bw.write("\t"+geneList);
//                }
//                else {
//                    bw.write("\t"+"-");
//                }
                //get gene 1000000 upstream and 5000 downstream (start-5000)(end+100000)
                // upstream (start, end+100000)
                List<Integer> geneIds = getGenesWithGeneCache(pos,pos+500000,chr);
                if (!geneIds.isEmpty()){
                    logger.debug("\t\tGetting genes Upstream 500000");
                    String geneList = listOfGenesToPrint(geneIds);
                    bw.write("\t"+geneList);
                }
                else {
                    bw.write("\t"+"-");
                }

                // downstream (start-5000, end) check if (start - 5000) < 0, set to 0
                int downstream = (pos<=500000) ? 0 : pos-500000;
                geneIds = getGenesWithGeneCache(downstream,pos,chr);
                if (!geneIds.isEmpty()){
                    logger.debug("\t\tGetting genes downstream 500000");
                    String geneList = listOfGenesToPrint(geneIds);
                    bw.write("\t"+geneList);
                }
                else {
                    bw.write("\t"+"-");
                }
                bw.write("\n");
            } // end while
            bw.close();
            br.close();
        }
        logger.info("Total runtime -- elapsed time: "+
                Utils.formatElapsedTime(pipeStart,System.currentTimeMillis()));
    }

    void run2() throws Exception{
        logger.info(version);
        long pipeStart = System.currentTimeMillis();
        logger.info("Pipeline started at "+sdt.format(new Date(pipeStart))+"\n");
        geneCacheMap = new HashMap<>();
        String file = "data/hrdp_118strains_plus_F1_genotype_all_chr.gvcf.gz";
        BufferedReader br = openFile(file);
        BufferedWriter bw = new BufferedWriter(new FileWriter("updated_hrdp_118strains_plus_F1_genotype_all_chr.gvcf"));
        String lineData;
        List<VariantMapData> tobeInserted = new ArrayList<>();
        List<VariantMapData> existing = new ArrayList<>();
        List<VariantSampleDetail> sampleDetails = new ArrayList<>();
        List<VariantMapData> updateGenicstatus= new ArrayList<>();
        while ((lineData = br.readLine()) != null) {
            if (lineData.startsWith("#")){
                bw.write(lineData+"\n");
                continue;
            }
            String[] cols = lineData.split("\t");
            // col[0] chr, col[1] pos, col[2] id (replace with rsId), col[3] ref, col[4] alt
            int pos = 0;
            String chr = "";
            String id="";
            String ref="";
            String alt="";
            String restOfLine = "";
            for (int i = 0; i < cols.length; i++){
                switch (i){
                    case 0:
                        chr = cols[i].replace("chr","");
                        break;
                    case 1:
                        pos = Integer.parseInt(cols[i]);
                        break;
                    case 2:
                        id = cols[i];
                        break;
                    case 3:
                        ref = cols[i].trim();
                        break;
                    case 4:
                        alt = cols[i].trim();
                        break;
                    default:
                        if (i+1== cols.length)
                            restOfLine += cols[i];
                        else
                            restOfLine += (cols[i]+"\t");
                        break;
                }
            } // end for
            if (Utils.isStringEmpty(chr) || pos == 0)
                continue;
            List<VariantMapData> vars = dao.getVariantsByLocation(chr,pos);
            boolean found = false;
            if (vars.isEmpty()){
                // create variant and replace id with rgdID
                VariantMapData v = createVariant(cols, tobeInserted);
                sampleDetails.add(createNewVariantSampleDetail(v));
                bw.write("chr"+chr+"\t"+pos+"\t"+v.getId()+"\t"+ref+"\t"+alt+"\t"+restOfLine);
            }
            else {
                if (vars.size() == 1) {
                    VariantMapData v = vars.get(0);
                    if (!Utils.isStringEmpty(v.getRsId()) && !Utils.stringsAreEqual(v.getRsId(), ".")) {
                        if (Utils.stringsAreEqual(v.getReferenceNucleotide(), ref) && Utils.stringsAreEqual(v.getVariantNucleotide(), alt)) {
                            id = v.getRsId();
                            found = true;
                            existing.add(v); // make samples
                        }
                    }
                    else {
                        id = v.getId() + "";
                    }
                    if (!found){
                        v = createVariant(cols,tobeInserted);
                    }
                    List<VariantSampleDetail> sd = dao.getVariantSampleDetail((int) v.getId(), sampleId);
                    if (sd.isEmpty()){
                        sampleDetails.add(createNewVariantSampleDetail(v));
                    }
                    String genicStatus = isGenic((int) v.getStartPos(), (int) v.getEndPos(), v.getChromosome()) ? "GENIC":"INTERGENIC";
                    if ( !Utils.stringsAreEqual(genicStatus, v.getGenicStatus()) || Utils.isStringEmpty(v.getGenicStatus()) ) {
                        v.setGenicStatus(genicStatus);
                        updateGenicstatus.add(v);
                    }
                    bw.write("chr" + chr + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + restOfLine);
                } else {
                    VariantMapData vmd = new VariantMapData();
                    for (VariantMapData variantMapData : vars) {
                        if (Utils.stringsAreEqual(variantMapData.getReferenceNucleotide(), ref) && Utils.stringsAreEqual(variantMapData.getVariantNucleotide(), alt)) {
                            vmd = variantMapData;
                            found=true;
                            existing.add(vmd); // make samples
                        }
                    }
                    if (!found){
                        vmd = createVariant(cols,tobeInserted);
                    }
                    if (vmd.getId() != 0 && !Utils.isStringEmpty(vmd.getRsId()) && !Utils.stringsAreEqual(vmd.getRsId(), ".")){
                        id = vmd.getRsId();
                    }
                    else
                        id = vmd.getId()+"";
                    List<VariantSampleDetail> sd = dao.getVariantSampleDetail((int) vmd.getId(), sampleId);
                    if (sd.isEmpty()){
                        sampleDetails.add(createNewVariantSampleDetail(vmd));
                    }

                    String genicStatus = isGenic((int) vmd.getStartPos(), (int) vmd.getEndPos(), vmd.getChromosome()) ? "GENIC":"INTERGENIC";
                    if ( !Utils.stringsAreEqual(genicStatus, vmd.getGenicStatus()) || Utils.isStringEmpty(vmd.getGenicStatus()) ) {
                        vmd.setGenicStatus(genicStatus);
                        updateGenicstatus.add(vmd);
                    }
                    bw.write("chr" + chr + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + restOfLine);

                }
            }
            bw.write("\n");
        }// end while

        br.close();
        bw.close();

        if (!updateGenicstatus.isEmpty()) {
            logger.info("\t\t\tVariants Genic Status being updated: "+updateGenicstatus.size());
            dao.updateGenicStatus(updateGenicstatus);
        }
        if (!tobeInserted.isEmpty()){
            logger.info("\tNew Variants being entered: " + tobeInserted.size());
            dao.insertVariants(tobeInserted);
            dao.insertVariantMapData(tobeInserted);
        }
        if (!sampleDetails.isEmpty()){
            logger.info("\tNew Sample Details being entered: "+ sampleDetails.size());
            dao.insertVariantSample(sampleDetails);
        }

        logger.info("Total runtime -- elapsed time: "+
                Utils.formatElapsedTime(pipeStart,System.currentTimeMillis()));
    }

    void closestGene() throws Exception {
        File folder = new File("data/human_snps/");
        List<File> files = listFilesInFolder(folder);
        geneCacheMap = new HashMap<>();
        logger.info(version);
        long pipeStart = System.currentTimeMillis();
        logger.info("Pipeline started at "+sdt.format(new Date(pipeStart))+"\n");
//        System.out.println(files.size());
        for (File file : files){
            BufferedReader br = openFile(file.getAbsolutePath());
            logger.info("\tGetting genes for: "+file.getName());
            BufferedWriter bw = new BufferedWriter(new FileWriter("closest_genes_from_"+file.getName()));
            String lineData;
            while ((lineData = br.readLine()) != null) {
                // get gene in snp region
                if (lineData.startsWith("Locus")) {
                    bw.write(lineData+"\tClose Gene\n");
                    continue;
                }
                bw.write(lineData);
                String[] cols = lineData.split("\t");
                String chr = cols[1];
                int pos = Integer.parseInt(cols[2]);
                int boundLimit = 0;
                List<Integer> geneIds = getClosestGene(boundLimit, pos, chr);
                while (geneIds.isEmpty()) {
                    boundLimit = boundLimit + 1000;
                    geneIds = getClosestGene(boundLimit, pos, chr);
                }

                String geneList = listOfGenesToPrint(geneIds);
                bw.write("\t"+geneList);

                bw.write("\n");
            } // end while
            bw.close();
            br.close();
        }
        logger.info("Total runtime -- elapsed time: "+
                Utils.formatElapsedTime(pipeStart,System.currentTimeMillis()));
    }

    List<File> listFilesInFolder(File folder) throws Exception {
        return Arrays.asList(Objects.requireNonNull(folder.listFiles()));
    }

    BufferedReader openFile(String fileName) throws IOException {

        String encoding = "UTF-8"; // default encoding

        InputStream is;
        if( fileName.endsWith(".gz") ) {
            is = new GZIPInputStream(new FileInputStream(fileName));
        } else {
            is = new FileInputStream(fileName);
        }

        BufferedReader reader = new BufferedReader(new InputStreamReader(is, encoding));
        return reader;
    }

    VariantMapData createVariant(String[] cols, List<VariantMapData> insert) throws Exception{
        String chr = "";
        int pos = 0;
        String id = "", ref = "", alt = "";
        for (int i = 0; i < cols.length; i++){
            switch (i){
                case 0:
                    chr = cols[i].replace("chr","");
                    break;
                case 1:
                    pos = Integer.parseInt(cols[i]);
                    break;
                case 2:
                    id = cols[i];
                    break;
                case 3:
                    ref = cols[i].trim();
                    break;
                case 4:
                    alt = cols[i].trim();
                    break;
                default:
                    break;
            }
        } // end for
        VariantMapData vmd = new VariantMapData();
        int speciesKey= SpeciesType.getSpeciesTypeKeyForMap(3);
        RgdId r = dao.createRgdId(RgdId.OBJECT_KEY_VARIANTS, "ACTIVE", "created by EVA pipeline", 372);
        vmd.setId(r.getRgdId());
        vmd.setSpeciesTypeKey(speciesKey);
        vmd.setVariantType("snp");
        vmd.setChromosome(chr);
        vmd.setStartPos(pos);
        vmd.setReferenceNucleotide(ref);
        vmd.setVariantNucleotide(alt);
        vmd.setEndPos(pos+1);
        vmd.setMapKey(372);
        String genicStat = isGenic((int) vmd.getStartPos(), (int) vmd.getEndPos(), vmd.getChromosome()) ? "GENIC":"INTERGENIC";
        vmd.setGenicStatus(genicStat);
        insert.add(vmd);
        return vmd;
    }

    List<Integer> getClosestGene(int boundLimit, int pos, String chr) throws Exception{
        GeneCache geneCache = geneCacheMap.get(chr);
        if( geneCache==null ) {
            geneCache = new GeneCache();
            geneCacheMap.put(chr, geneCache);
            geneCache.loadCache(38, chr, DataSourceFactory.getInstance().getDataSource());
        }
        return geneCache.calculateClosestGene(boundLimit,pos);
    }


    boolean isGenic(int start, int stop, String chr) throws Exception {

        GeneCache geneCache = geneCacheMap.get(chr);
        if( geneCache==null ) {
            geneCache = new GeneCache();
            geneCacheMap.put(chr, geneCache);
            geneCache.loadCache(38, chr, DataSourceFactory.getInstance().getDataSource());
        }
        List<Integer> geneRgdIds = geneCache.getGeneRgdIds(start,stop);
        return !geneRgdIds.isEmpty();
    }
    List<Integer> getGenesWithGeneCache(int start, int stop, String chr) throws Exception {

       GeneCache geneCache = geneCacheMap.get(chr);
        if( geneCache==null ) {
            geneCache = new GeneCache();
            geneCacheMap.put(chr, geneCache);
            geneCache.loadCache(38, chr, DataSourceFactory.getInstance().getDataSource());
        }
        List<Integer> geneRgdIds = geneCache.getGeneRgdIds(start,stop);
        return geneRgdIds;
    }

    public VariantSampleDetail createNewVariantSampleDetail(VariantMapData vmd) throws Exception{
        VariantSampleDetail vsd = new VariantSampleDetail(); // add to variant_sample_detail with eva sample leave zygosity stuff empty
        vsd.setId(vmd.getId());
        vsd.setSampleId(sampleId);
        vsd.setDepth(9);
        vsd.setVariantFrequency(1);
        return vsd;
    }

    String listOfGenesToPrint(List<Integer> geneIds) throws Exception{
        StringBuilder genes = new StringBuilder();
        List<String> geneList = dao.getGeneNamesInRegion(geneIds);
        if (!geneList.isEmpty()){
            for (int i = 0; i < geneList.size(); i++){
                if (i == (geneList.size()-1))
                    genes.append(geneList.get(i));
                else
                    genes.append(geneList.get(i)).append(", ");
            }
            return genes.toString();
        }
        else
            return "-";

    }

    Map<String, GeneCache> geneCacheMap;

    public void setVersion(String version) {
        this.version = version;
    }

    public String getVersion() {
        return version;
    }

    public void setSampleId(Integer sampleId) {
        this.sampleId = sampleId;
    }

    public Integer getSampleId() {
        return sampleId;
    }
}