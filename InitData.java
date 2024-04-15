package EconomicDispatch;


import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;
import java.util.HashMap;


public class InitData {
    public static void main(String[] args) {
        InitData initData = new InitData();
        initData.init();
    }

    public HashMap<String,double[]> init(){
        String originDataFile = "system\\AllData.txt";
        int N_TGU = 13,N_CHP = 6,N_HOU = 3,N_SPVP = 3,N_WTG = 3;
        int Dim = N_TGU + N_CHP*2 + N_HOU + N_SPVP + N_WTG;
        double[] a_TGU = new double[N_TGU], b_TGU = new double[N_TGU], c_TGU = new double[N_TGU],d_TGU = new double[N_TGU], e_TGU = new double[N_TGU], UR_TGU = new double[N_TGU], DR_TGU = new double[N_TGU];
        double[] a_CHP = new double[N_CHP], b_CHP = new double[N_CHP], c_CHP = new double[N_CHP],d_CHP = new double[N_CHP], e_CHP = new double[N_CHP], f_CHP = new double[N_CHP], UR_CHP = new double[N_CHP], DR_CHP = new double[N_CHP];
        double[] a_HOU = new double[N_HOU], b_HOU = new double[N_HOU], c_HOU = new double[N_HOU];
        double[] S_SPVP = new double[1+24], T_ref = new double[1+24], T_amb = new double[1+24], P_SPVP_Rated = new double[N_SPVP], a_SPVP = new double[N_SPVP], CC_SPVP = new double[N_SPVP];
        double[] V_Wind_Real = new double[1+24],V_Wind_CUT_IN = new double[1+24],V_Wind_CUT_OFF=new double[1+24],V_Wind_Rated = new double[1+24], CC_WTG = new double[N_WTG],P_Wind_Rated = new double[N_WTG];
        double[] PD = new double[1+24], HD = new double[1+24];
        int TGU_SIZE = 0,CHP_SIZE = 0,HOU_SIZE = 0, t = 1;
        double[] P_CHP_MAX = new double[N_CHP], P_CHP_MIN = new double[N_CHP], P_TGU_MAX = new double[N_TGU], P_TGU_MIN = new double[N_TGU];
        double[] H_CHP_MAX = new double[N_CHP], H_CHP_MIN = new double[N_CHP], H_HOU_MAX = new double[N_HOU], H_HOU_MIN = new double[N_HOU];
        double[] H_CHP_MID = new double[N_CHP];
        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(originDataFile));
            String s = br.readLine();
            while(s != null){
                if(s.equals("TGU")){
                    s = br.readLine();
                    while(s != null && !s.equals("CHP")){
                        String[] data = s.split(" ");
                        P_TGU_MIN[TGU_SIZE] = Double.parseDouble(data[1]);
                        P_TGU_MAX[TGU_SIZE] = Double.parseDouble(data[2]);
                        a_TGU[TGU_SIZE] = Double.parseDouble(data[3]);
                        b_TGU[TGU_SIZE] = Double.parseDouble(data[4]);
                        c_TGU[TGU_SIZE] = Double.parseDouble(data[5]);
                        d_TGU[TGU_SIZE] = Double.parseDouble(data[6]);
                        e_TGU[TGU_SIZE] = Double.parseDouble(data[7]);
                        UR_TGU[TGU_SIZE] = Double.parseDouble(data[8]);
                        DR_TGU[TGU_SIZE] = Double.parseDouble(data[9]);
                        TGU_SIZE++;
                        s = br.readLine();
                    }
                }else if(s.equals("CHP")){
                    s = br.readLine();
                    while(s != null && !s.equals("HOU")){
                        String[] data = s.split(" ");
                        P_CHP_MIN[CHP_SIZE] = Double.parseDouble(data[1]);
                        P_CHP_MAX[CHP_SIZE] = Double.parseDouble(data[2]);
                        H_CHP_MIN[CHP_SIZE] = Double.parseDouble(data[3]);
                        H_CHP_MAX[CHP_SIZE] = Double.parseDouble(data[4]);
                        H_CHP_MID[CHP_SIZE] = Double.parseDouble(data[5]);
                        a_CHP[CHP_SIZE] = Double.parseDouble(data[6]);
                        b_CHP[CHP_SIZE] = Double.parseDouble(data[7]);
                        c_CHP[CHP_SIZE] = Double.parseDouble(data[8]);
                        d_CHP[CHP_SIZE] = Double.parseDouble(data[9]);
                        e_CHP[CHP_SIZE] = Double.parseDouble(data[10]);
                        f_CHP[CHP_SIZE] = Double.parseDouble(data[11]);
                        UR_CHP[CHP_SIZE] = Double.parseDouble(data[12]);
                        DR_CHP[CHP_SIZE] = Double.parseDouble(data[13]);
                        CHP_SIZE++;
                        s = br.readLine();
                    }
                }else if(s.equals("HOU")){
                    s = br.readLine();
                    while(s != null && !s.equals("Solar")){
                        String[] data = s.split(" ");
                        H_HOU_MIN[HOU_SIZE] = Double.parseDouble(data[1]);
                        H_HOU_MAX[HOU_SIZE] = Double.parseDouble(data[2]);
                        a_HOU[HOU_SIZE] = Double.parseDouble(data[3]);
                        b_HOU[HOU_SIZE] = Double.parseDouble(data[4]);
                        c_HOU[HOU_SIZE] = Double.parseDouble(data[5]);
                        HOU_SIZE++;
                        s = br.readLine();
                    }
                } else if (s.equals("Solar")) {
                    s = br.readLine();
                    while(s != null && !s.equals("Wind")){
                        String[] data = s.split(" ");
                        S_SPVP[t] = Double.parseDouble(data[1]);
                        T_ref[t] = Double.parseDouble(data[2]);
                        T_amb[t] = 24;
                        t++;
                        s = br.readLine();
                    }
                    for(int x = 0; x < N_SPVP; x++){
                        P_SPVP_Rated[x] = 175;
                        a_SPVP[x] = 0.005;
                        CC_SPVP[x] = 7;
                    }
                } else if (s.equals("Wind")) {
                    s = br.readLine();
                    t = 1;
                    while(s != null && !s.equals("Need")){
                        String[] data = s.split(" ");
                        V_Wind_Real[t] = Double.parseDouble(data[2])*0.5 + Double.parseDouble(data[1])*0.5;
                        s = br.readLine();
                        t++;
                    }
                    for(int x = 0; x < N_WTG; x++){
                        V_Wind_Rated[x] = 15;
                        V_Wind_CUT_OFF[x] = 25;
                        V_Wind_CUT_IN[x] = 4;
                        P_Wind_Rated[x] = 175;
                        CC_WTG[x] = 6;
                    }
                }else {
                    s = br.readLine();
                    t = 1;
                    while(s != null){
                        String[] data = s.split(" ");
                        PD[t] = Double.parseDouble(data[1]);
                        HD[t] = Double.parseDouble(data[2]);
                        t++;
                        s = br.readLine();
                    }
                }
            }
            br.close();
        }catch (Exception ignored){

        }
        HashMap<String, double[]> ret = new HashMap<>();
        double[] basicInfo = new double[5];
        basicInfo[0] = N_TGU;basicInfo[1] = N_CHP;basicInfo[2] = N_HOU;basicInfo[3] = N_SPVP;basicInfo[4] = N_WTG;
        ret.put("basicInfo",basicInfo);ret.put("P_Need",PD);ret.put("H_Need",HD);
        ret.put("a_TGU",a_TGU);ret.put("b_TGU",b_TGU);ret.put("c_TGU",c_TGU);ret.put("d_TGU",d_TGU);ret.put("e_TGU",e_TGU);ret.put("UR_TGU",UR_TGU);ret.put("DR_TGU",DR_TGU);
        ret.put("a_CHP",a_CHP);ret.put("b_CHP",b_CHP);ret.put("c_CHP",c_CHP);ret.put("d_CHP",d_CHP);ret.put("e_CHP",e_CHP);ret.put("f_CHP",f_CHP);ret.put("UR_CHP",UR_CHP);ret.put("DR_CHP",DR_CHP);
        ret.put("a_HOU",a_HOU);ret.put("b_HOU",b_HOU);ret.put("c_HOU",c_HOU);
        ret.put("T_ref",T_ref);ret.put("T_amb",T_amb);ret.put("S_SPVP",S_SPVP);ret.put("a_SPVP",a_SPVP);ret.put("P_SPVP_Rated",P_SPVP_Rated);
        ret.put("V_Wind_Real",V_Wind_Real);ret.put("V_Wind_CUT_IN",V_Wind_CUT_IN);ret.put("V_Wind_CUT_OFF",V_Wind_CUT_OFF);ret.put("V_Wind_Rated",V_Wind_Rated);ret.put("P_Wind_Rated",P_Wind_Rated);

        ret.put("P_CHP_MAX",P_CHP_MAX);ret.put("P_CHP_MIN",P_CHP_MIN);ret.put("P_TGU_MIN",P_TGU_MIN);ret.put("P_TGU_MAX",P_TGU_MAX);
        ret.put("H_CHP_MIN",H_CHP_MIN);ret.put("H_CHP_MAX",H_CHP_MAX);ret.put("H_HOU_MIN",H_HOU_MIN);ret.put("H_HOU_MAX",H_HOU_MAX);

        ret.put("CC_SPVP",CC_SPVP);ret.put("CC_WTG",CC_WTG);
        ret.put("H_CHP_MID",H_CHP_MID);
        return ret;
    }
}
