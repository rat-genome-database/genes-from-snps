package edu.mcw.rgd.genesFromSnps;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
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

                }
            }
        }
        catch (Exception e) {
            Utils.printStackTrace(e, LogManager.getLogger("status"));
            throw e;
        }
    }

    void run() throws Exception {
        File folder = new File("data/");
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
                    bw.write(lineData+"\tGene\tGenes Upstream\tGenes Downstream\n");
                    continue;
                }
                bw.write(lineData);
                String[] cols = lineData.split("\t");
                String chr = cols[1];
                int pos = Integer.parseInt(cols[2]);
                List<Integer> geneIds = getGenesWithGeneCache(pos,pos,chr);
                if (!geneIds.isEmpty()){
                    logger.debug("\t\tGetting Gene");
                    String geneList = listOfGenesToPrint(geneIds);
                    bw.write("\t"+geneList);
                }
                else {
                    bw.write("\t"+"-");
                }
                //get gene 1000000 upstream and 5000 downstream (start-5000)(end+100000)
                // upstream (start, end+100000)
                geneIds = getGenesWithGeneCache(pos,pos+100000,chr);
                if (!geneIds.isEmpty()){
                    logger.debug("\t\tGetting genes Upstream 100000");
                    String geneList = listOfGenesToPrint(geneIds);
                    bw.write("\t"+geneList);
                }
                else {
                    bw.write("\t"+"-");
                }

                // downstream (start-5000, end) check if (start - 5000) < 0, set to 0
                int downstream = (pos<=5000) ? 0 : pos-5000;
                geneIds = getGenesWithGeneCache(downstream,pos,chr);
                if (!geneIds.isEmpty()){
                    logger.debug("\t\tGetting genes downstream 5000");
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

        String file = "data/hrdp_118strains_plus_F1_genotype_all_chr.gvcf.gz";
        BufferedReader br = openFile(file);
        BufferedWriter bw = new BufferedWriter(new FileWriter("updated_hrdp_118strains_plus_F1_genotype_all_chr.gvcf"));
        String lineData;

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
                        ref = cols[i];
                        break;
                    case 4:
                        alt = cols[i];
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
            if (vars.isEmpty()){
                bw.write("chr"+chr+"\t"+pos+"\t"+id+"\t"+ref+"\t"+alt+"\t"+restOfLine);
            }
            else {
                if (vars.size() == 1) {
                    VariantMapData v = vars.get(0);
                    if (!Utils.isStringEmpty(v.getRsId()) && !Utils.stringsAreEqual(v.getRsId(), ".")) {
                        if (Utils.stringsAreEqual(v.getReferenceNucleotide(), ref) && Utils.stringsAreEqual(v.getVariantNucleotide(), alt)) {
                            id = v.getRsId();
                        }
                    }
                    bw.write("chr" + chr + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + restOfLine);
                } else {
                    VariantMapData vmd = new VariantMapData();
                    for (VariantMapData variantMapData : vars) {
                        if (Utils.stringsAreEqual(variantMapData.getReferenceNucleotide(), ref) && Utils.stringsAreEqual(variantMapData.getVariantNucleotide(), alt)) {
                            vmd = variantMapData;
                        }
                    }
                    if (vmd.getId() != 0 && !Utils.isStringEmpty(vmd.getRsId()) && !Utils.stringsAreEqual(vmd.getRsId(), "."))
                        id = vmd.getRsId();

                    bw.write("chr" + chr + "\t" + pos + "\t" + id + "\t" + ref + "\t" + alt + "\t" + restOfLine);

                }
            }
            bw.write("\n");
        }// end while

        br.close();
        bw.close();
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
}