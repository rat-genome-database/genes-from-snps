package edu.mcw.rgd.genesFromSnps;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.*;
import edu.mcw.rgd.dao.impl.variants.VariantDAO;
import edu.mcw.rgd.dao.spring.variants.VariantMapQuery;
import edu.mcw.rgd.dao.spring.variants.VariantSampleQuery;
import edu.mcw.rgd.datamodel.*;
import edu.mcw.rgd.datamodel.ontologyx.Term;
import edu.mcw.rgd.datamodel.variants.VariantMapData;
import edu.mcw.rgd.datamodel.variants.VariantSampleDetail;
import edu.mcw.rgd.process.Utils;
import oracle.jdbc.proxy.annotation.Pre;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.springframework.jdbc.core.SqlParameter;
import org.springframework.jdbc.object.BatchSqlUpdate;

import javax.sql.DataSource;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.Types;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Created by llamers on 1/28/2020.
 */
public class DAO {

    private OntologyXDAO xdao = new OntologyXDAO();
    private GeneDAO geneDAO = new GeneDAO();
    private MapDAO mdao = new MapDAO();
    private RGDManagementDAO managementDAO = new RGDManagementDAO();
    private VariantDAO vdao = new VariantDAO();

    public String getConnection(){
        return mdao.getConnectionInfo();
    }

    public DataSource getVariantDataSource() throws Exception{
        return DataSourceFactory.getInstance().getCarpeNovoDataSource();
    }

    List<String> getGeneNamesInRegion(List<Integer> geneRgdIds) throws Exception{
        List<String> genes = new ArrayList<>();
        String geneIds = "";
        for (int i = 0; i < geneRgdIds.size(); i++){
            if ( (geneRgdIds.size()-1) == i){
                geneIds += geneRgdIds.get(i);
            }
            else {
                geneIds += geneRgdIds.get(i)+", ";
            }

        }
        String sql = "SELECT g.GENE_SYMBOL FROM GENES g, RGD_IDS r WHERE r.object_status='ACTIVE' AND r.rgd_id=g.rgd_id AND r.RGD_ID in ("+geneIds+")";
        Connection c = null;
        try{
            c = geneDAO.getConnection();
            PreparedStatement ps = c.prepareStatement(sql);
            ResultSet rs = ps.executeQuery();
            while (rs.next()){
                genes.add(rs.getString(1));
            }
        }
        catch (Exception e){
            e.printStackTrace();
        }
        finally {
            c.close();
        }
        return genes;
    }

    List<Gene> getGenesInRegion(int start, int stop, String chr) throws Exception{
        List<MapData> mapData = mdao.getMapDataWithinRange(start,stop,chr,38,0);
        List<Gene> geneList = new ArrayList<>();
        if (mapData.size()>0) {
            GeneDAO gdao = new GeneDAO();
            for (MapData m : mapData) {
                Gene g = gdao.getGene(m.getRgdId());
                if (g != null)
                    geneList.add(g);
            }
        }
        return geneList;
    }

    List<VariantMapData> getVariantsByLocation(String chr, int pos) throws Exception{
        return vdao.getVariantsWithGeneLocation(372,chr,pos,pos);
    }

    public RgdId createRgdId(int objectKey, String objectStatus, String notes, int mapKey) throws Exception{
        int speciesKey=SpeciesType.getSpeciesTypeKeyForMap(mapKey);
        return managementDAO.createRgdId(objectKey, objectStatus, notes, speciesKey);
    }

    public void insertVariants(List<VariantMapData> mapsData)  throws Exception{
        BatchSqlUpdate sql1 = new BatchSqlUpdate(this.getVariantDataSource(),
                "INSERT INTO variant (RGD_ID,REF_NUC, VARIANT_TYPE, VAR_NUC, RS_ID, CLINVAR_ID, SPECIES_TYPE_KEY)" +
                        "VALUES (?,?,?,?,?,?,?)",
                new int[]{Types.INTEGER,Types.VARCHAR,Types.VARCHAR, Types.VARCHAR, Types.VARCHAR, Types.VARCHAR,Types.INTEGER}, 10000);
        sql1.compile();
        for( VariantMapData v: mapsData) {
            long id = v.getId();
            sql1.update(id, v.getReferenceNucleotide(), v.getVariantType(), v.getVariantNucleotide(), v.getRsId(), v.getClinvarId(), v.getSpeciesTypeKey());

        }
        sql1.flush();
    }
    public void insertVariantMapData(List<VariantMapData> mapsData)  throws Exception{
        BatchSqlUpdate sql2 = new BatchSqlUpdate(this.getVariantDataSource(),
                "INSERT INTO variant_map_data (RGD_ID,CHROMOSOME,START_POS,END_POS,PADDING_BASE,GENIC_STATUS,MAP_KEY) " +
                        "VALUES ( ?,?,?,?,?,?,?)",
                new int[]{Types.INTEGER,Types.VARCHAR, Types.INTEGER, Types.INTEGER, Types.VARCHAR,Types.VARCHAR, Types.INTEGER});
        sql2.compile();
        for( VariantMapData v: mapsData) {
            long id = v.getId();
            sql2.update(id, v.getChromosome(), v.getStartPos(), v.getEndPos(), v.getPaddingBase(), v.getGenicStatus(), v.getMapKey());
        }
        sql2.flush();
    }

    public List<VariantSampleDetail> getVariantSampleDetail(int rgdId, int sampleId) throws Exception{
        String sql = "SELECT * FROM variant_sample_detail  WHERE rgd_id=? AND sample_id=?";
        VariantSampleQuery q = new VariantSampleQuery(getVariantDataSource(), sql);
        q.declareParameter(new SqlParameter(Types.INTEGER));
        q.declareParameter(new SqlParameter(Types.INTEGER));
        return q.execute(rgdId, sampleId);
    }

    public void insertVariantSample(List<VariantSampleDetail> sampleData) throws Exception {
        BatchSqlUpdate bsu= new BatchSqlUpdate(this.getVariantDataSource(),
                "INSERT INTO variant_sample_detail (RGD_ID,SAMPLE_ID,TOTAL_DEPTH,VAR_FREQ) " +
                        "VALUES (?,?,?,?)",
                new int[]{Types.INTEGER,Types.INTEGER, Types.INTEGER, Types.INTEGER});
        bsu.compile();
        for(VariantSampleDetail v: sampleData ) {
            bsu.update(v.getId(), v.getSampleId(),v.getDepth(),v.getVariantFrequency());
        }
        bsu.flush();
    }

    public void insertSampleDetail(VariantSampleDetail vs) throws Exception{
        vdao.insertSampleDetail(vs);
    }

    public void updateGenicStatus(List<VariantMapData> mapsData) throws Exception {
        BatchSqlUpdate sql2 = new BatchSqlUpdate(this.getVariantDataSource(),
                "update variant_map_data set GENIC_STATUS=? where RGD_ID=?",
                new int[]{Types.VARCHAR,Types.INTEGER});
        sql2.compile();
        for( VariantMapData v: mapsData) {
            long id = v.getId();
            sql2.update(v.getGenicStatus(),id);
        }
        sql2.flush();
    }

}
