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
            Main main =(Main) bf.getBean("main");
            main.run();
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
                    logger.info("\t\tGetting Gene");
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
                    logger.info("\t\tGetting genes Upstream 100000");
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
                    logger.info("\t\tGetting genes downstream 5000");
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