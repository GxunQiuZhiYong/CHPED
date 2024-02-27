package EconomicDispatch;

import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;

public class SCAALG {

    public static class SCA{
        public static Random random = new Random();
        public static int Dim,popSize,T,N_TGU,N_CHP,N_HOU,N_SPVP,N_WTG;
        public static double[] a_TGU,a_CHP,a_HOU,b_TGU,b_CHP,b_HOU,c_TGU,c_CHP,c_HOU,d_TGU,e_TGU,d_CHP,e_CHP,f_CHP,CC_SPVP,CC_WTG,UR_TGU,DR_TGU,UR_CHP,DR_CHP;
        public static double[][][] P_H_MAX,P_H_MIN;
        // 太阳能模型
        public static double[] T_ref, T_amb, S_SPVP, a_SPVP, P_SPVP_Reted;
        // 风电模型
        public static double[] V_Wind_Real,V_Wind_CUT_IN,V_Wind_CUT_OFF,V_Wind_Rated,P_Wind_Rated;
        // 约束
        public static double[] P_Need,H_Need;
        public static double[] P_CHP_MIN, P_CHP_MAX;
        public static double[][][] p,prep;
        public static double[] fitVal;
        public static double[][] pbest;
        public static double pbestFitness;
        public static double windCost,SpvpCost;
        public static double[] H_CHP_MID;
        public SCA(int popSize, int T){
            SCA.popSize = popSize;
            SCA.T = T;
        }
        public static void windPower(){
            windCost = 0;
            for(int t = 1; t <= 24; t++){
                for(int i = 0;i < N_WTG; i++){
                    double tmp = 0;
                    if(V_Wind_Real[t] < V_Wind_CUT_IN[i] || V_Wind_Real[t] > V_Wind_CUT_OFF[i]) tmp = 0;
                    else if (V_Wind_Real[t] >= V_Wind_CUT_IN[i] && V_Wind_Real[t] <= V_Wind_Rated[i]) {
                        tmp = P_Wind_Rated[i] * (V_Wind_Real[t] - V_Wind_CUT_IN[i])/(V_Wind_Rated[i] - V_Wind_CUT_IN[i]);
                    }else {
                        tmp = P_Wind_Rated[i];
                    }
                    windCost += CC_WTG[i]*tmp;
                    P_Need[t] -= tmp;
                }
            }
        }
        public static void solarPower(){
            SpvpCost = 0;
            for(int t = 1; t <= 24; t++){
                for(int i = 0; i < N_SPVP; i++){
                    double tmp = Math.min(Math.max(P_SPVP_Reted[i]*(1+a_SPVP[i]*(T_ref[t]-T_amb[t]))*S_SPVP[t]/1000,0),P_SPVP_Reted[i]);
                    P_Need[t] -= tmp;
                    SpvpCost += tmp*CC_SPVP[i];
                }
            }
        }
        public static boolean constrainP(int s,int t){
            double p_total = 0;
            for(int i = 0; i < Dim - N_HOU;i++){
                p_total += p[s][t][i];
                if(i >= N_TGU) i++;
            }
            if(P_Need[t]*1.02 > p_total) {
                return false;
            }
            return true;
        }
        public static boolean constrainH(int s, int t){
            double h_total = 0;
            for(int i = N_TGU;i < Dim;i++){
                if(i < N_TGU + N_CHP*2){
                    i++;
                    h_total += p[s][t][i];
                }else h_total += p[s][t][i];
            }
            if(H_Need[t]*1.02 > h_total) {
                return false;
            }
            return true;
        }
        public static boolean constrain(int s){
            for(int t = 1; t <= 24; t++){
                double p_total = 0, h_total = 0;
                for(int i = 0; i < Dim; i++){
                    if(i >= Dim - N_HOU || (i >= N_TGU && i < N_TGU + N_CHP*2 && (i - N_TGU)%2 == 1)) h_total += p[s][t][i];
                    else p_total += p[s][t][i];
                }
                if(p_total - P_Need[t]*1.02 < 0 || h_total - H_Need[t]*1.02 < 0) return false;
            }
            return true;
        }
        public static void CHP_FOR(int s){
            for(int t = 1; t <= 24; t++){
                int b = 0;
                for(int i = N_TGU; i < N_CHP*2 + N_TGU; i+=2){
                    P_H_MIN[s][t][i] = Math.max(0.75*(p[s][t][i+1] - H_CHP_MID[b])+(P_CHP_MIN[b]-0.15*H_CHP_MID[b]), P_CHP_MIN[b] - 0.15*p[s][t][i+1]);
                    P_H_MAX[s][t][i] = Math.max(P_CHP_MAX[b] - 0.15*p[s][t][i+1],0);
                    b++;
                }
            }
        }
        public static void updateH(int s, int t){
            int b = 0;
            for(int i = N_TGU; i < Dim; i++){
                if(i < N_TGU + N_CHP * 2) i++;
                prep[s][t][i] = p[s][t][i];
                if(i < Dim - N_HOU){
                    P_H_MIN[s][t][i-1] = Math.max(0.75*(p[s][t][i] - H_CHP_MID[b])+(P_CHP_MIN[b]-0.15*H_CHP_MID[b]), P_CHP_MIN[b] - 0.15*p[s][t][i]);
                    P_H_MAX[s][t][i-1] = P_CHP_MAX[b] - 0.15*p[s][t][i];
                }
                b++;
            }
        }
        public static void updateP(int s){
            for(int t = 1; t <= 24; t++){
                for(int i = 0; i < Dim-N_HOU; i++){
                    prep[s][t][i] = p[s][t][i];
                    if(i >= N_TGU && i < N_TGU + N_CHP * 2) i++;
                }
            }

        }
        public static void backPreH(int s, int t){
            int b = 0;
            for(int i = N_TGU; i < Dim; i++){
                if(i < N_TGU + N_CHP * 2) i++;
                p[s][t][i] = prep[s][t][i];
                if(i < Dim - N_HOU){
                    P_H_MIN[s][t][i-1] = Math.max(0.75*(p[s][t][i] - H_CHP_MID[b])+(P_CHP_MIN[b]-0.15*H_CHP_MID[b]), P_CHP_MIN[b] - 0.15*p[s][t][i]);
                    P_H_MAX[s][t][i-1] = P_CHP_MAX[b] - 0.15*p[s][t][i];
                }
                b++;
            }
        }
        public static void backPreP(int s){
            for(int t = 1;t <= 24; t++){
                for(int i = 0; i < Dim - N_HOU; i++){
                    p[s][t][i] = prep[s][t][i];
                    if(i >= N_TGU && i < N_TGU + N_CHP * 2) i++;
                }
            }
        }
        public static double getFitVal(int s){
            double ret = 0;
            int a = 0, b = 0,c = 0;
            for(int i = 0,j = 0; i < Dim; i++,j++){
                for(int t = 1; t <= 24; t++){
                    if(i < N_TGU){
                        ret += a_TGU[a] + b_TGU[a]*p[s][t][i] + c_TGU[a]*Math.pow(p[s][t][i],2) + Math.abs(d_TGU[a]*Math.sin(e_TGU[a]*(p[s][t][i] - P_H_MIN[s][t][i])));
                    } else if (i >= N_TGU && i < N_TGU + N_CHP*2) {
                        ret += a_CHP[b] + b_CHP[b]*p[s][t][i] + c_CHP[b]*Math.pow(p[s][t][i],2) + d_CHP[b]*p[s][t][i+1] + e_CHP[b]*Math.pow(p[s][t][i+1],2) + f_CHP[b]*p[s][t][i+1]*p[s][t][i];
                    } else if (i >= N_CHP*2 + N_TGU) {
                        ret += a_HOU[c] + b_HOU[c]*p[s][t][i] + c_HOU[c]*Math.pow(p[s][t][i],2);
                    }
                }

                if(i < N_TGU) a++;
                else if (i >= N_TGU && i < N_TGU + N_CHP*2){
                    b++;
                    i++;
                }
                else c++;
            }
            return ret;
        }
        public static void fitness(){
            for(int s = 0; s < popSize; s++){
                double fi = getFitVal(s);
                if(pbestFitness > fi){
                    pbestFitness = fi;
                    for(int x = 1; x <= 24; x++){
                        for(int y = 0; y < N_TGU + N_CHP*2 + N_HOU; y++){
                            pbest[x][y] = p[s][x][y];
                        }
                    }
                }
                fitVal[s] = fi;
            }
        }
        public static void updatePStage(double r1, double r2, double r3, double r4){
            for(int s = 0; s < popSize; s++){
                // 更新H
                for(int t = 1; t <= 24; t++){
                    for(int i = N_TGU; i < Dim; i++){
                        if(i < N_TGU + N_CHP * 2) i++;
                        if(r4 > 0.5) p[s][t][i] = Math.min(Math.max(p[s][t][i]+r1*Math.sin(r2)*Math.abs(r3*pbest[t][i]-p[s][t][i]), P_H_MIN[s][t][i]), P_H_MAX[s][t][i]);
                        else p[s][t][i] = Math.min(Math.max(p[s][t][i]+r1*Math.cos(r2)*Math.abs(r3*pbest[t][i]-p[s][t][i]), P_H_MIN[s][t][i]), P_H_MAX[s][t][i]);
                    }
                    if(constrainH(s,t)) updateH(s,t);
                    else backPreH(s,t);
                }

                // 更新P
                boolean flag = true;
                for(int t = 1; t <= 24; t++){
                    if(!flag) continue;
                    int a = 0,bb = 0;
                    for(int i = 0; i < Dim - N_HOU; i++){
                        if(r4 > 0.5) p[s][t][i] = Math.min(Math.max(p[s][t][i]+r1*Math.sin(r2)*Math.abs(r3*pbest[t][i]-p[s][t][i]), P_H_MIN[s][t][i]), P_H_MAX[s][t][i]);
                        else p[s][t][i] = Math.min(Math.max(p[s][t][i]+r1*Math.cos(r2)*Math.abs(r3*pbest[t][i]-p[s][t][i]), P_H_MIN[s][t][i]), P_H_MAX[s][t][i]);
                        if(t >= 2){
                            if(i < N_TGU){
                                if(p[s][t][i] - p[s][t-1][i] > 0) p[s][t][i] = Math.min(p[s][t][i], p[s][t-1][i] + UR_TGU[a]);
                                else p[s][t][i] = Math.max(p[s][t][i],p[s][t-1][i]-DR_TGU[a]);
                                a++;
                            } else if (i < N_CHP*2 + N_TGU && (i - N_TGU)%2==0) {
                                if(p[s][t][i] - p[s][t-1][i] > 0) p[s][t][i] = Math.min(p[s][t][i], p[s][t-1][i] + UR_CHP[bb]);
                                else p[s][t][i] = Math.max(p[s][t][i],p[s][t-1][i]-DR_CHP[bb]);
                                bb++;
                            }
                        }
                        if(i >= N_TGU) i++;
                    }
                    if(!constrainP(s,t)){
                        backPreP(s);
                        flag = false;
                    }
                }
                if(flag) updateP(s);

            }
        }
        public static void init() {
            InitData initData = new InitData();
            Map<String,double[]> data = initData.init();
            double[] basicInfo = data.get("basicInfo");
            N_TGU = (int) basicInfo[0];N_CHP = (int) basicInfo[1];N_HOU = (int) basicInfo[2];N_SPVP = (int) basicInfo[3];N_WTG = (int) basicInfo[4];
            Dim = N_TGU + N_CHP*2 + N_HOU;
            a_TGU = data.get("a_TGU");b_TGU = data.get("b_TGU");c_TGU = data.get("c_TGU");d_TGU = data.get("d_TGU");e_TGU = data.get("e_TGU");
            a_CHP = data.get("a_CHP");b_CHP = data.get("b_CHP");c_CHP = data.get("c_CHP");d_CHP = data.get("d_CHP");e_CHP = data.get("e_CHP");f_CHP = data.get("f_CHP");
            a_HOU = data.get("a_HOU");b_HOU = data.get("b_HOU");c_HOU = data.get("c_HOU");
            CC_SPVP = data.get("CC_SPVP");CC_WTG = data.get("CC_WTG");
            UR_TGU = data.get("UR_TGU");DR_TGU = data.get("DR_TGU");UR_CHP = data.get("UR_CHP");DR_CHP = data.get("DR_CHP");
            P_CHP_MAX = data.get("P_CHP_MAX");P_CHP_MIN = data.get("P_CHP_MIN");
            P_H_MAX = new double[popSize][1+24][Dim];
            P_H_MIN = new double[popSize][1+24][Dim];
            p = new double[popSize][1+24][Dim];prep = new double[popSize][1+24][Dim];
            fitVal = new double[popSize];
            pbest = new double[1+24][N_TGU + N_CHP*2 + N_HOU];
            pbestFitness = 1e40;
            for(int s = 0; s < popSize; s++){
                for(int t = 1; t <= 24; t++){
                    P_H_MAX[s][t] = data.get("P_H_MAX");
                    P_H_MIN[s][t] = data.get("P_H_MIN");
                }
                fitVal[s] = 1e40;
            }
            T_ref = data.get("T_ref");T_amb = data.get("T_amb");a_SPVP = data.get("a_SPVP");P_SPVP_Reted = data.get("P_SPVP_Rated");S_SPVP = data.get("S_SPVP");
            V_Wind_CUT_IN = data.get("V_Wind_CUT_IN");V_Wind_CUT_OFF = data.get("V_Wind_CUT_OFF");V_Wind_Rated = data.get("V_Wind_Rated");P_Wind_Rated = data.get("P_Wind_Rated");
            V_Wind_Real = data.get("V_Wind_Real");
            P_Need = data.get("P_Need");H_Need = data.get("H_Need");
            H_CHP_MID = data.get("H_CHP_MID");
            windPower();
            solarPower();
            System.out.println(Arrays.toString(P_Need));
            System.out.println(Arrays.toString(H_Need));
            for(int s = 0; s < popSize; s++){
                System.out.println(s);
                do{
                    for(int i = 0,j = 0; i < Dim; i++,j++){
                        if(i >= N_TGU && i < N_TGU + 2*N_CHP) i++;
                        for(int t = 1; t <= 24; t++){
                            if(t >= 2 && i < N_TGU){
                                p[s][t][i] = random.nextDouble(Math.max(p[s][t-1][i] - DR_TGU[i],P_H_MIN[s][t][i]), Math.min(p[s][t-1][i] + UR_TGU[i],P_H_MAX[s][t][i]));
                            }
                            else p[s][t][i] = random.nextDouble(P_H_MAX[s][t][i]*0.5,P_H_MAX[s][t][i]);
                        }
                    }
                    CHP_FOR(s);
                    for(int i = N_TGU,j = 0; i < N_TGU + 2*N_CHP;i+=2,j++){
                        for(int t = 1; t <= 24; t++){
                            if(t >= 2){
                                p[s][t][i] = random.nextDouble(Math.max(p[s][t-1][i] - DR_CHP[j],P_H_MIN[s][t][i]), Math.min(p[s][t-1][i] + UR_CHP[j],P_H_MAX[s][t][i]));
                            }
                            else p[s][t][i] = random.nextDouble(P_H_MIN[s][t][i],P_H_MAX[s][t][i]);
                        }
                    }
                }while (!constrain(s));
                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < Dim; i++){
                        prep[s][t][i] = p[s][t][i];
                    }
                }
            }
            fitness();
        }
        public static void myQuicksort(int low, int high){
            if(low < high){
                int l = low, r = high;
                double tmp = fitVal[l];
                double[][] tmp1 = p[l];
                while(l < r){
                    while(l < r && fitVal[r] > tmp) r--;
                    if(l < r){
                        fitVal[l] = fitVal[r];
                        p[l] = p[r];
                        l++;
                    }
                    while(l < r && fitVal[l] < tmp) l++;
                    if(l < r){
                        fitVal[r] = fitVal[l];
                        p[r] = p[l];
                        r--;
                    }
                }
                fitVal[l] = tmp;
                p[l] = tmp1;
                myQuicksort(low, l-1);
                myQuicksort(l+1, high);
            }
        }
        public void start(String filename) {
            double a = 2;
            FileOutputStream foss = null;
            try {
                foss = new FileOutputStream("C:\\Users\\Qiu\\Desktop\\EconomicDispatch\\result3\\system4\\bestits.txt",true);
            }catch (Exception e){

            }
            init();
            for(int its = 1; its <= T; its++){
                myQuicksort(0,popSize-1);
                updatePStage(a - a*its/T,random.nextDouble()*2*Math.PI, random.nextDouble(0,10), random.nextDouble());
                fitness();
                if(its % 1 == 0) System.out.println("iterations : " + its + " / " + T + "  | minFit : " + (pbestFitness + windCost + SpvpCost));
                try {
                    foss.write((String.valueOf(pbestFitness+ windCost + SpvpCost) + '\n').getBytes());

                }catch (Exception e){
                }
            }
            try {
                foss.close();
                FileOutputStream fos = new FileOutputStream(filename,true);
                String s = String.valueOf(pbestFitness+ windCost + SpvpCost) + " |  " + Arrays.deepToString(pbest) + '\n';
                fos.write(s.getBytes());
                fos.close();
            }catch (Exception ignored){

            }
        }
    }

    public static void main(String[] args) {
        for(int eachtime = 0 ; eachtime < 30; eachtime++){
            SCA sca = new SCA(100,1000);
            String filename = "C:\\Users\\Qiu\\Desktop\\EconomicDispatch\\result3\\system4\\data.txt";
            sca.start(filename);
        }
    }
}


