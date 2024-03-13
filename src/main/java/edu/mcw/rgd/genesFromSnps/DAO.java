package edu.mcw.rgd.genesFromSnps;

import edu.mcw.rgd.dao.DataSourceFactory;
import edu.mcw.rgd.dao.impl.*;
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

    public String getConnection(){
        return mdao.getConnectionInfo();
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

}
