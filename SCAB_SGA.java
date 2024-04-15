package EconomicDispatch;

import java.io.FileOutputStream;
import java.util.Arrays;
import java.util.Map;
import java.util.Random;

public class SCAB_SGA {

    public static class SGA{
        public static Random random = new Random();
        public static int popSize,T,N_TGU,N_CHP,N_HOU,N_SPVP,N_WTG;
        public static double[] a_TGU,a_CHP,a_HOU,b_TGU,b_CHP,b_HOU,c_TGU,c_CHP,c_HOU,d_TGU,e_TGU,d_CHP,e_CHP,f_CHP,CC_SPVP,CC_WTG,UR_TGU,DR_TGU,UR_CHP,DR_CHP;
        // 太阳能模型
        public static double[] T_ref, T_amb, S_SPVP, a_SPVP, P_SPVP_Reted;
        // 风电模型
        public static double[] V_Wind_Real,V_Wind_CUT_IN,V_Wind_CUT_OFF,V_Wind_Rated,P_Wind_Rated;
        // 约束
        public static double[] P_Need,H_Need,fitVal;
        public static double[] P_CHP_MIN, P_CHP_MAX, H_CHP_MIN, H_CHP_MAX, P_TGU_MIN, P_TGU_MAX, H_HOU_MAX, H_HOU_MIN, P_P2G_MAX,H_CHP_MID;

        public static double[][][] P_CHP_CUR_MAX,P_CHP_CUR_MIN;
        public static double[][][] p_TGU,p_CHP,h_CHP,h_HOU,p_SPVP,p_WTG,p_P2G, prep_TGU,prep_CHP,preh_CHP,preh_HOU,prep_SPVP,prep_WTG,prep_P2G,v_TGU,v_p_CHP,v_h_CHP,v_SPVP,v_WTG,v_HOU,v_P2G;
        public static double[][] P_SPVP_MAX, P_WTG_MAX;
        public static double[][] best_TGU, best_p_CHP, best_h_CHP, best_HOU, best_SPVP, best_WTG, best_P2G;
        public static double pbestFitness;
        public static double r1,r3,r4;
        public SGA(int popSize, int T, double r1, double r3, double r4){
            SGA.popSize = popSize;
            SGA.T = T;
            SGA.r1 = r1;
            SGA.r3 = r3;
            SGA.r4 = r4;
        }

        public static void windPower(){
            for(int t = 1; t <= 24; t++){
                for(int i = 0;i < N_WTG; i++){
                    double tmp;
                    if(V_Wind_Real[t] < V_Wind_CUT_IN[i] || V_Wind_Real[t] > V_Wind_CUT_OFF[i]) tmp = 0;
                    else if (V_Wind_Real[t] >= V_Wind_CUT_IN[i] && V_Wind_Real[t] <= V_Wind_Rated[i]) {
                        tmp = P_Wind_Rated[i] * (V_Wind_Real[t] - V_Wind_CUT_IN[i])/(V_Wind_Rated[i] - V_Wind_CUT_IN[i]);
                    }else {
                        tmp = P_Wind_Rated[i];
                    }
                    P_WTG_MAX[t][i] = tmp;
                }
            }
        }
        public static void solarPower(){
            for(int t = 1; t <= 24; t++){
                for(int i = 0; i < N_SPVP; i++){
                    double tmp = Math.min(Math.max(P_SPVP_Reted[i]*(1+a_SPVP[i]*(T_ref[t]-T_amb[t]))*S_SPVP[t]/1000,0),P_SPVP_Reted[i]);
                    P_SPVP_MAX[t][i] = tmp;
                }
            }
        }
        public static double SPVP_punish(int s, int t, int i){
            return 2*GetIntegrals.get("SPVP_Oe",0,p_SPVP[s][t][i]) + GetIntegrals.get("SPVP_Ue", p_SPVP[s][t][i], P_SPVP_MAX[t][i]);

        }
        public static double WTG_punish(int s, int t, int i){
            return 2*GetIntegrals.get("WTG_Oe",0,p_WTG[s][t][i]) + GetIntegrals.get("WTG_Ue", p_WTG[s][t][i], P_WTG_MAX[t][i]);
        }
        public static void CHP_FOR(int s){
            for(int t = 1; t <= 24; t++){
                for(int i = 0; i < N_CHP; i++){
                    P_CHP_CUR_MIN[s][t][i] = Math.max(0.75*(h_CHP[s][t][i] - H_CHP_MID[i])+(P_CHP_MIN[i]-0.15*H_CHP_MID[i]), P_CHP_MIN[i] - 0.15*h_CHP[s][t][i])-p_P2G[s][t][i]*(1+1.02*0.5);
                    P_CHP_CUR_MAX[s][t][i] = P_CHP_MAX[i] - 0.15*h_CHP[s][t][i];
                }
            }
        }
        public static double getFitVal(int s){
            double ret = 0.0;
            for(int t = 1; t <= 24; t++){
                for(int i = 0; i < N_TGU; i++){
                    ret += a_TGU[i] + b_TGU[i]*p_TGU[s][t][i] + c_TGU[i]*Math.pow(p_TGU[s][t][i],2) + Math.abs(d_TGU[i]*Math.sin(e_TGU[i]*(p_TGU[s][t][i] - P_TGU_MIN[i])));
                }
                for(int i = 0; i < N_CHP; i++){
                    ret += a_CHP[i] + b_CHP[i]*p_CHP[s][t][i] + c_CHP[i]*Math.pow(p_CHP[s][t][i],2) + d_CHP[i]*h_CHP[s][t][i] + e_CHP[i]*Math.pow(h_CHP[s][t][i],2) + f_CHP[i]*p_CHP[s][t][i]*h_CHP[s][t][i];
                }
                for(int i = 0; i < N_HOU; i++){
                    ret += a_HOU[i] + b_HOU[i]*h_HOU[s][t][i] + c_HOU[i]*Math.pow(h_HOU[s][t][i],2);
                }
                for(int i = 0; i < N_SPVP; i++){
                    ret += CC_SPVP[i]*p_SPVP[s][t][i] + SPVP_punish(s,t,i);
                }
                for(int i = 0; i < N_WTG; i++){
                    ret += CC_WTG[i]*p_WTG[s][t][i] + WTG_punish(s,t,i);
                }
            }
            return ret;
        }
        public static void fitness(){
            for(int s = 0; s < popSize; s++){
                double f = getFitVal(s);
                if(pbestFitness > f){
                    pbestFitness = f;
                    for(int t = 1; t <= 24; t++){
                        System.arraycopy(p_TGU[s][t], 0, best_TGU[t], 0, N_TGU);
                        System.arraycopy(p_CHP[s][t], 0, best_p_CHP[t], 0, N_CHP);
                        System.arraycopy(h_HOU[s][t], 0, best_HOU[t], 0, N_HOU);
                        System.arraycopy(p_SPVP[s][t], 0, best_SPVP[t], 0, N_SPVP);
                        System.arraycopy(p_WTG[s][t], 0, best_WTG[t], 0, N_WTG);
                        System.arraycopy(h_CHP[s][t], 0, best_h_CHP[t], 0, N_CHP);
                        System.arraycopy(p_P2G[s][t], 0, best_P2G[t], 0, N_CHP);
                    }
                }
                fitVal[s] = f;
            }
        }
        public static boolean constrainP(int s,int t){
            double p_total = 0;
            for(int i = 0; i < N_TGU;i++) p_total += p_TGU[s][t][i];
            for(int i = 0; i < N_CHP;i++) p_total += (p_CHP[s][t][i] + p_P2G[s][t][i]*0.4*0.55);
            for(int i = 0; i < N_SPVP;i++) p_total += p_SPVP[s][t][i];
            for(int i = 0; i < N_WTG;i++) p_total += p_WTG[s][t][i];
            return !(P_Need[t] * 1.1 > p_total);
        }
        public static boolean constrainH(int s, int t){
            double h_total = 0;
            for(int i = 0;i < N_CHP;i++) h_total += h_CHP[s][t][i];
            for(int i = 0;i < N_HOU;i++) h_total += h_HOU[s][t][i];
            return !(H_Need[t] * 1.15 > h_total);
        }
        public static void updateH(int s, int t){
            for(int i = 0; i < N_CHP; i++){
                preh_CHP[s][t][i] = h_CHP[s][t][i];
                P_CHP_CUR_MIN[s][t][i] = Math.max(0.75*(h_CHP[s][t][i] - H_CHP_MID[i])+(P_CHP_MIN[i]-0.15*H_CHP_MID[i]), P_CHP_MIN[i] - 0.15*h_CHP[s][t][i])-p_P2G[s][t][i]*(1+1.02*0.5);
                P_CHP_CUR_MAX[s][t][i] = P_CHP_MAX[i] - 0.15*h_CHP[s][t][i];
            }
            System.arraycopy(h_HOU[s][t], 0, preh_HOU[s][t], 0, N_HOU);
        }
        public static void updateP(int s){
            for(int t = 1; t <= 24; t++){
                System.arraycopy(p_TGU[s][t], 0, prep_TGU[s][t], 0, N_TGU);
                System.arraycopy(p_CHP[s][t], 0, prep_CHP[s][t], 0, N_CHP);
                System.arraycopy(p_SPVP[s][t], 0, prep_SPVP[s][t], 0, N_SPVP);
                System.arraycopy(p_WTG[s][t], 0, prep_WTG[s][t], 0, N_WTG);
                System.arraycopy(p_P2G[s][t], 0, prep_P2G[s][t], 0, N_CHP);
            }
        }
        public static void backPreH(int s, int t){
            for(int i = 0; i < N_CHP; i++){
                h_CHP[s][t][i] = preh_CHP[s][t][i];
                P_CHP_CUR_MIN[s][t][i] = Math.max(0.75*(h_CHP[s][t][i] - H_CHP_MID[i])+(P_CHP_MIN[i]-0.15*H_CHP_MID[i]), P_CHP_MIN[i] - 0.15*h_CHP[s][t][i])-p_P2G[s][t][i]*(1+1.02*0.5);
                P_CHP_CUR_MAX[s][t][i] = P_CHP_MAX[i] - 0.15*h_CHP[s][t][i];
            }
            System.arraycopy(preh_HOU[s][t], 0, h_HOU[s][t], 0, N_HOU);
        }
        public static void backPreP(int s){
            for(int t = 1;t <= 24; t++){
                System.arraycopy(prep_TGU[s][t], 0, p_TGU[s][t], 0, N_TGU);
                System.arraycopy(prep_CHP[s][t], 0, p_CHP[s][t], 0, N_CHP);
                System.arraycopy(prep_SPVP[s][t], 0, p_SPVP[s][t], 0, N_SPVP);
                System.arraycopy(prep_WTG[s][t], 0, p_WTG[s][t], 0, N_WTG);
                System.arraycopy(prep_P2G[s][t], 0, p_P2G[s][t], 0, N_CHP);
            }
        }
        public static void updateV(double its, double xita){
            for(int s = 0; s < popSize; s++){
                for(int t = 1; t <= 24; t++){
                    // H
                    for(int i = 0; i < N_CHP; i++){
                        v_h_CHP[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_h_CHP[s][t][i] + (best_h_CHP[t][i] - h_CHP[s][t][i] - 1.29 * v_h_CHP[s][t][i] * v_h_CHP[s][t][i] * Math.sin(xita)) * 0.01;
                        if(Math.abs(v_h_CHP[s][t][i]) > H_CHP_MAX[i]) {
                            v_h_CHP[s][t][i] = random.nextDouble()*H_CHP_MAX[i]*(v_h_CHP[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_h_CHP[s][t][i]) < H_CHP_MAX[i]/T && random.nextDouble() < 0.8) {
                            v_h_CHP[s][t][i] *= 10;
                        }
                    }
                    for(int i = 0; i < N_HOU; i++){
                        v_HOU[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_HOU[s][t][i] + (best_HOU[t][i] - h_HOU[s][t][i] - 1.29 * v_HOU[s][t][i] * v_HOU[s][t][i] * Math.sin(xita)) * 0.01;
                        if(Math.abs(v_HOU[s][t][i]) > H_HOU_MAX[i]) {
                            v_HOU[s][t][i] = random.nextDouble()*H_HOU_MAX[i]*(v_HOU[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_HOU[s][t][i]) < H_HOU_MAX[i]/T && random.nextDouble() < 0.8) {
                            v_HOU[s][t][i] *= 10;
                        }
                    }
                    // P
                    for(int i = 0; i < N_TGU; i++){
                        v_TGU[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_TGU[s][t][i] + (best_TGU[t][i] - p_TGU[s][t][i] - 1.29 * v_TGU[s][t][i] * v_TGU[s][t][i] * Math.sin(xita)) * 0.01;
                        if(Math.abs(v_TGU[s][t][i]) > P_TGU_MAX[i]) {
                            v_TGU[s][t][i] = random.nextDouble()*P_TGU_MAX[i]*(v_TGU[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_TGU[s][t][i]) < P_TGU_MAX[i]/T*r3 && random.nextDouble() < r1) {
                            v_TGU[s][t][i] *= 10;
                        }
                    }
                    for(int i = 0; i < N_CHP; i++){
                        v_p_CHP[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_p_CHP[s][t][i] + (best_p_CHP[t][i] - p_CHP[s][t][i] - 1.29 * v_p_CHP[s][t][i] * v_p_CHP[s][t][i] * Math.sin(xita)) * 0.01;
                        v_P2G[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_P2G[s][t][i] + (best_P2G[t][i] - p_P2G[s][t][i] - 1.29 * v_P2G[s][t][i] * v_P2G[s][t][i] * Math.sin(xita)) * 0.01;
                        if(Math.abs(v_p_CHP[s][t][i]) > P_CHP_CUR_MAX[s][t][i]) {
                            v_p_CHP[s][t][i] = random.nextDouble()*P_CHP_CUR_MAX[s][t][i]*(v_p_CHP[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_p_CHP[s][t][i]) < P_CHP_CUR_MAX[s][t][i]/T*r3 && random.nextDouble() < r1) {
                            v_p_CHP[s][t][i] *= 10;
                        }
                        if(Math.abs(v_P2G[s][t][i]) > P_P2G_MAX[i]) {
                            v_P2G[s][t][i] = random.nextDouble()*P_P2G_MAX[i]*(v_P2G[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_P2G[s][t][i]) < P_P2G_MAX[i]/T*r3 && random.nextDouble() < r1) {
                            v_P2G[s][t][i] *= 10;
                        }
                    }
                    for(int i = 0; i < N_SPVP; i++){
                        v_SPVP[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_SPVP[s][t][i] + (best_SPVP[t][i] - p_SPVP[s][t][i] - 1.29 * v_SPVP[s][t][i] * v_SPVP[s][t][i] * Math.sin(xita)) * 0.01;
                        if(Math.abs(v_SPVP[s][t][i]) > P_SPVP_MAX[t][i]) {
                            v_SPVP[s][t][i] = random.nextDouble()*P_SPVP_MAX[t][i]*(v_SPVP[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_SPVP[s][t][i]) < P_SPVP_MAX[t][i]/T*r3 && random.nextDouble() < r1) {
                            v_SPVP[s][t][i] *= 10;
                        }
                    }
                    for(int i = 0; i < N_WTG; i++){
                        v_WTG[s][t][i] = 4 * (its / T) / Math.pow(Math.E, (4 * its / T)) * v_WTG[s][t][i] + (best_WTG[t][i] - p_WTG[s][t][i] - 1.29 * v_WTG[s][t][i] * v_WTG[s][t][i] * Math.sin(xita)) * 0.01;
                        if(Math.abs(v_WTG[s][t][i]) > P_WTG_MAX[t][i]) {
                            v_WTG[s][t][i] = random.nextDouble()*P_WTG_MAX[t][i]*(v_WTG[s][t][i] > 0 ? 1 : -1);
                        }
                        else if (Math.abs(v_WTG[s][t][i]) < P_WTG_MAX[t][i]/T*r3 && random.nextDouble() < r1) {
                            v_WTG[s][t][i] *= 10;
                        }
                    }
                }
            }
        }
        public static void updateStageOne(double b, double d, double e){
            double[][] pc_TGU = new double[25][N_TGU], pc_p_CHP = new double[25][N_CHP], pc_h_CHP = new double[25][N_CHP], pc_HOU = new double[25][N_HOU],pc_SPVP = new double[25][N_SPVP], pc_WTG = new double[25][N_WTG], pc_P2G = new double[25][N_CHP];
            double sum = Arrays.stream(fitVal).sum()*popSize;
            for(int s = 0; s < popSize; s++){
                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < N_TGU; i++){
                        pc_TGU[t][i] += fitVal[s]*p_TGU[s][t][i];
                        if(s == popSize -1) pc_TGU[t][i] = pc_TGU[t][i]/sum;
                    }
                    for(int i = 0; i < N_CHP; i++){
                        pc_p_CHP[t][i] += fitVal[s]*p_CHP[s][t][i];
                        if(s == popSize -1) pc_p_CHP[t][i] = pc_p_CHP[t][i]/sum;
                    }
                    for(int i = 0; i < N_CHP; i++){
                        pc_h_CHP[t][i] += fitVal[s]*h_CHP[s][t][i];
                        pc_P2G[t][i] += fitVal[s]*p_P2G[s][t][i];
                        if(s == popSize -1) {
                            pc_h_CHP[t][i] = pc_h_CHP[t][i]/sum;
                            pc_P2G[t][i] = pc_P2G[t][i]/sum;
                        }
                    }
                    for(int i = 0; i < N_HOU; i++){
                        pc_HOU[t][i] += fitVal[s]*h_HOU[s][t][i];
                        if(s == popSize -1) pc_HOU[t][i] = pc_HOU[t][i]/sum;
                    }
                    for(int i = 0; i < N_SPVP; i++){
                        pc_SPVP[t][i] += fitVal[s]*p_SPVP[s][t][i];
                        if(s == popSize -1) pc_SPVP[t][i] = pc_SPVP[t][i]/sum;
                    }
                    for(int i = 0; i < N_WTG; i++){
                        pc_WTG[t][i] += fitVal[s]*p_WTG[s][t][i];
                        if(s == popSize -1) pc_WTG[t][i] = pc_WTG[t][i]/sum;
                    }
                }
            }
            for(int s = 0; s < popSize; s++){
                // H
                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < N_CHP; i++){
                        if(s < popSize*0.2){
                            h_CHP[s][t][i] = Math.min(Math.max(h_CHP[s][t][i] + b*(best_h_CHP[t][i]-h_CHP[s][t][i]) + v_h_CHP[s][t][i], H_CHP_MIN[i]), H_CHP_MAX[i]);
                            p_P2G[s][t][i] = Math.min(Math.max(p_P2G[s][t][i] + b*(best_P2G[t][i]-p_P2G[s][t][i]) + v_P2G[s][t][i], 0), P_P2G_MAX[i]);
                        } else if (s > popSize*0.8) {
                            h_CHP[s][t][i] = Math.min(Math.max(h_CHP[s][t][i] + b*(best_h_CHP[t][i]-h_CHP[s][t][i]) + d*(pc_h_CHP[t][i]-h_CHP[s][t][i]) + v_h_CHP[s][t][i], H_CHP_MIN[i]), H_CHP_MAX[i]);
                            p_P2G[s][t][i] = Math.min(Math.max(p_P2G[s][t][i] + b*(best_P2G[t][i]-p_P2G[s][t][i]) + d*(pc_P2G[t][i]-p_P2G[s][t][i]) + v_P2G[s][t][i], 0), P_P2G_MAX[i]);
                        }else {
                            h_CHP[s][t][i] = Math.min(Math.max(h_CHP[s][t][i] + b*(best_h_CHP[t][i]-h_CHP[s][t][i]) + d*(pc_h_CHP[t][i]-h_CHP[s][t][i]) - e*(h_CHP[popSize-1][t][i]+h_CHP[s][t][i]) + v_h_CHP[s][t][i], H_CHP_MIN[i]), H_CHP_MAX[i]);
                            p_P2G[s][t][i] = Math.min(Math.max(p_P2G[s][t][i] + b*(best_P2G[t][i]-p_P2G[s][t][i]) + d*(pc_P2G[t][i]-p_P2G[s][t][i]) + v_P2G[s][t][i], 0), P_P2G_MAX[i]);
                        }
                    }
                    for(int i = 0; i < N_HOU; i++){
                        if(s < popSize*0.2){
                            h_HOU[s][t][i] = Math.min(Math.max(h_HOU[s][t][i] + b*(best_HOU[t][i]-h_HOU[s][t][i]) + v_HOU[s][t][i], H_HOU_MIN[i]), H_HOU_MAX[i]);
                        } else if (s > popSize*0.8) {
                            h_HOU[s][t][i] = Math.min(Math.max(h_HOU[s][t][i] + b*(best_HOU[t][i]-h_HOU[s][t][i]) + d*(pc_HOU[t][i]-h_HOU[s][t][i]) + v_HOU[s][t][i], H_HOU_MIN[i]), H_HOU_MAX[i]);
                        }else {
                            h_HOU[s][t][i] = Math.min(Math.max(h_HOU[s][t][i] + b*(best_HOU[t][i]-h_HOU[s][t][i]) + d*(pc_HOU[t][i]-h_HOU[s][t][i]) - e*(h_HOU[popSize-1][t][i]+h_HOU[s][t][i]) + v_HOU[s][t][i], H_HOU_MIN[i]), H_HOU_MAX[i]);
                        }
                    }
                    if(constrainH(s,t)) updateH(s,t);
                    else backPreH(s,t);
                }

                // P
                boolean flag = true;
                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < N_TGU; i++){
                        if(s < popSize*0.2){
                            p_TGU[s][t][i] = Math.min(Math.max(p_TGU[s][t][i] + b*(best_TGU[t][i]-p_TGU[s][t][i]) + v_TGU[s][t][i], P_TGU_MIN[i]), P_TGU_MAX[i]);
                        } else if (s > popSize*0.8) {
                            p_TGU[s][t][i] = Math.min(Math.max(p_TGU[s][t][i] + b*(best_TGU[t][i]-p_TGU[s][t][i]) + d*(pc_TGU[t][i]-p_TGU[s][t][i]) + v_TGU[s][t][i], P_TGU_MIN[i]), P_TGU_MAX[i]);
                        }else {
                            p_TGU[s][t][i] = Math.min(Math.max(p_TGU[s][t][i] + b*(best_TGU[t][i]-p_TGU[s][t][i]) + d*(pc_TGU[t][i]-p_TGU[s][t][i]) - e*(p_TGU[popSize-1][t][i]+p_TGU[s][t][i]) + v_TGU[s][t][i], P_TGU_MIN[i]), P_TGU_MAX[i]);
                        }

                        if(t >= 2){
                            if(p_TGU[s][t][i] - p_TGU[s][t-1][i] > 0) p_TGU[s][t][i] = Math.min(p_TGU[s][t][i], p_TGU[s][t-1][i] + UR_TGU[i]);
                            else p_TGU[s][t][i] = Math.max(p_TGU[s][t][i],p_TGU[s][t-1][i]-DR_TGU[i]);

                            if(p_TGU[s][t][i] - p_TGU[s][t-1][i] > 0) p_TGU[s][t][i] = Math.min(p_TGU[s][t][i], p_TGU[s][t-1][i] + UR_TGU[i]);
                            else p_TGU[s][t][i] = Math.max(p_TGU[s][t][i],p_TGU[s][t-1][i]-DR_TGU[i]);
                        }
                    }
                    for(int i = 0; i < N_CHP; i++){
                        if(s < popSize*0.2){
                            p_CHP[s][t][i] = Math.min(Math.max(p_CHP[s][t][i] + b*(best_p_CHP[t][i]-p_CHP[s][t][i]) + v_p_CHP[s][t][i], P_CHP_CUR_MIN[s][t][i]), P_CHP_CUR_MAX[s][t][i]);
                        } else if (s > popSize*0.8) {
                            p_CHP[s][t][i] = Math.min(Math.max(p_CHP[s][t][i] + b*(best_p_CHP[t][i]-p_CHP[s][t][i]) + d*(pc_P2G[t][i]-p_CHP[s][t][i]) + v_p_CHP[s][t][i], P_CHP_CUR_MIN[s][t][i]), P_CHP_CUR_MAX[s][t][i]);
                        }else {
                            p_CHP[s][t][i] = Math.min(Math.max(p_CHP[s][t][i] + b*(best_p_CHP[t][i]-p_CHP[s][t][i]) + d*(pc_P2G[t][i]-p_CHP[s][t][i]) - e*(p_CHP[popSize-1][t][i]+p_CHP[s][t][i]) + v_p_CHP[s][t][i], P_CHP_CUR_MIN[s][t][i]), P_CHP_CUR_MAX[s][t][i]);
                        }
                        if(t >= 2){
                            if(p_CHP[s][t][i] - p_CHP[s][t-1][i] > 0) p_CHP[s][t][i] = Math.min(p_CHP[s][t][i], p_CHP[s][t-1][i] + UR_CHP[i]);
                            else p_CHP[s][t][i] = Math.max(p_CHP[s][t][i],p_CHP[s][t-1][i]-DR_CHP[i]);

                            if(p_CHP[s][t][i] - p_CHP[s][t-1][i] > 0) p_CHP[s][t][i] = Math.min(p_CHP[s][t][i], p_CHP[s][t-1][i] + UR_CHP[i]);
                            else p_CHP[s][t][i] = Math.max(p_CHP[s][t][i],p_CHP[s][t-1][i]-DR_CHP[i]);
                        }

                    }
                    for(int i = 0; i < N_SPVP; i++){
                        if(s < popSize*0.2){
                            p_SPVP[s][t][i] = Math.min(Math.max(p_SPVP[s][t][i] + b*(best_SPVP[t][i]-p_SPVP[s][t][i]) + v_SPVP[s][t][i], 0), P_SPVP_MAX[t][i]);
                        } else if (s > popSize*0.8) {
                            p_SPVP[s][t][i] = Math.min(Math.max(p_SPVP[s][t][i] + b*(best_SPVP[t][i]-p_SPVP[s][t][i]) + d*(pc_SPVP[t][i]-p_SPVP[s][t][i]) + v_SPVP[s][t][i], 0), P_SPVP_MAX[t][i]);
                        }else {
                            p_SPVP[s][t][i] = Math.min(Math.max(p_SPVP[s][t][i] + b*(best_SPVP[t][i]-p_SPVP[s][t][i]) + d*(pc_SPVP[t][i]-p_SPVP[s][t][i]) - e*(p_SPVP[popSize-1][t][i]+p_SPVP[s][t][i]) + v_SPVP[s][t][i], 0), P_SPVP_MAX[t][i]);
                        }
                    }
                    for(int i = 0; i < N_WTG; i++){
                        if(s < popSize*0.2){
                            p_WTG[s][t][i] = Math.min(Math.max(p_WTG[s][t][i] + b*(best_WTG[t][i]-p_WTG[s][t][i]) + v_WTG[s][t][i], 0), P_WTG_MAX[t][i]);
                        } else if (s > popSize*0.8) {
                            p_WTG[s][t][i] = Math.min(Math.max(p_WTG[s][t][i] + b*(best_WTG[t][i]-p_WTG[s][t][i]) + d*(pc_WTG[t][i]-p_WTG[s][t][i]) + v_WTG[s][t][i], 0), P_WTG_MAX[t][i]);
                        }else {
                            p_WTG[s][t][i] = Math.min(Math.max(p_WTG[s][t][i] + b*(best_WTG[t][i]-p_WTG[s][t][i]) + d*(pc_WTG[t][i]-p_WTG[s][t][i]) - e*(p_WTG[popSize-1][t][i]+p_WTG[s][t][i]) + v_WTG[s][t][i], 0), P_WTG_MAX[t][i]);
                        }
                    }
                    if(!constrainP(s,t)){
                        backPreP(s);
                        flag = false;
                    }
                }
                if(flag) updateP(s);
            }
        }
        public static double brown(double its){
            double d = 0.0;
            for(int i = 0; i < N_TGU + N_CHP*3 + N_HOU + N_SPVP + N_WTG; i++) d+= random.nextGaussian(0,2*(Math.cos((its%100) / 100 * (Math.PI/2)) + r4));
            return Math.sqrt((1.0d/(N_TGU + N_CHP*3 + N_HOU + N_SPVP + N_WTG)))*d;
        }
        public static void updateStageTwo(double its){
            double o = brown(its);
            for(int s = 0; s < popSize; s++){
                // H
                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < N_CHP; i++){
                        double r = random.nextDouble(0,1);
                        if(r > 0.5){
                            h_CHP[s][t][i] = Math.min(Math.max(h_CHP[s][t][i] + (h_CHP[s][t][i] - best_h_CHP[t][i])*r, H_CHP_MIN[i]), H_CHP_MAX[i]);
                            p_P2G[s][t][i] = Math.min(Math.max(p_P2G[s][t][i] + (p_P2G[s][t][i] - best_P2G[t][i])*r, 0), P_P2G_MAX[i]);
                        }else {
                            h_CHP[s][t][i] = Math.min(Math.max(best_h_CHP[t][i] + (h_CHP[s][t][i]-best_h_CHP[t][i])*r*o, H_CHP_MIN[i]), H_CHP_MAX[i]);
                            p_P2G[s][t][i] = Math.min(Math.max(best_P2G[t][i] + (p_P2G[s][t][i]-best_P2G[t][i])*r*o, 0), P_P2G_MAX[i]);
                        }
                    }
                    for(int i = 0; i < N_HOU; i++){
                        double r = random.nextDouble(0,1);
                        if(r > 0.5){
                            h_HOU[s][t][i] = Math.min(Math.max(h_HOU[s][t][i] + (h_HOU[s][t][i] - best_HOU[t][i])*r, H_HOU_MIN[i]), H_HOU_MAX[i]);
                        }else {
                            h_HOU[s][t][i] = Math.min(Math.max(best_HOU[t][i] + (h_HOU[s][t][i]-best_HOU[t][i])*r*o, H_HOU_MIN[i]), H_HOU_MAX[i]);
                        }
                    }
                    if(constrainH(s,t)) updateH(s,t);
                    else backPreH(s,t);
                }
                // P
                boolean flag = true;
                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < N_TGU; i++){
                        double r = random.nextDouble(0,1);
                        if(r > 0.5){
                            p_TGU[s][t][i] = Math.min(Math.max(p_TGU[s][t][i] + (p_TGU[s][t][i] - best_TGU[t][i])*r, P_TGU_MIN[i]), P_TGU_MAX[i]);
                        }else {
                            p_TGU[s][t][i] = Math.min(Math.max(best_TGU[t][i] + (p_TGU[s][t][i]-best_TGU[t][i])*r*o, P_TGU_MIN[i]), P_TGU_MAX[i]);
                        }
                        if(t >= 2){
                            if(p_TGU[s][t][i] - p_TGU[s][t-1][i] > 0) p_TGU[s][t][i] = Math.min(p_TGU[s][t][i], p_TGU[s][t-1][i] + UR_TGU[i]);
                            else p_TGU[s][t][i] = Math.max(p_TGU[s][t][i],p_TGU[s][t-1][i]-DR_TGU[i]);

                            if(p_TGU[s][t][i] - p_TGU[s][t-1][i] > 0) p_TGU[s][t][i] = Math.min(p_TGU[s][t][i], p_TGU[s][t-1][i] + UR_TGU[i]);
                            else p_TGU[s][t][i] = Math.max(p_TGU[s][t][i],p_TGU[s][t-1][i]-DR_TGU[i]);
                        }
                    }
                    for(int i = 0; i < N_CHP; i++){
                        double r = random.nextDouble(0,1);
                        if(r > 0.5){
                            p_CHP[s][t][i] = Math.min(Math.max(p_CHP[s][t][i] + (p_CHP[s][t][i] - best_p_CHP[t][i])*r, P_CHP_CUR_MIN[s][t][i]), P_CHP_CUR_MAX[s][t][i]);
                        }else {
                            p_CHP[s][t][i] = Math.min(Math.max(best_p_CHP[t][i] + (p_CHP[s][t][i]-best_p_CHP[t][i])*r*o, P_CHP_CUR_MIN[s][t][i]), P_CHP_CUR_MAX[s][t][i]);
                        }

                        if(t >= 2){
                            if(p_CHP[s][t][i] - p_CHP[s][t-1][i] > 0) p_CHP[s][t][i] = Math.min(p_CHP[s][t][i], p_CHP[s][t-1][i] + UR_CHP[i]);
                            else p_CHP[s][t][i] = Math.max(p_CHP[s][t][i],p_CHP[s][t-1][i]-DR_CHP[i]);

                            if(p_CHP[s][t][i] - p_CHP[s][t-1][i] > 0) p_CHP[s][t][i] = Math.min(p_CHP[s][t][i], p_CHP[s][t-1][i] + UR_CHP[i]);
                            else p_CHP[s][t][i] = Math.max(p_CHP[s][t][i],p_CHP[s][t-1][i]-DR_CHP[i]);
                        }
                    }
                    for(int i = 0; i < N_SPVP; i++){
                        double r = random.nextDouble(0,1);
                        if(r > 0.5){
                            p_SPVP[s][t][i] = Math.min(Math.max(p_SPVP[s][t][i] + (p_SPVP[s][t][i] - best_SPVP[t][i])*r, 0), P_SPVP_MAX[t][i]);
                        }else {
                            p_SPVP[s][t][i] = Math.min(Math.max(best_SPVP[t][i] + (p_SPVP[s][t][i]-best_SPVP[t][i])*r*o, 0), P_SPVP_MAX[t][i]);
                        }
                    }
                    for(int i = 0; i < N_WTG; i++){
                        double r = random.nextDouble(0,1);
                        if(r > 0.5){
                            p_WTG[s][t][i] = Math.min(Math.max(p_WTG[s][t][i] + (p_WTG[s][t][i] - best_WTG[t][i])*r, 0), P_WTG_MAX[t][i]);
                        }else {
                            p_WTG[s][t][i] = Math.min(Math.max(best_WTG[t][i] + (p_WTG[s][t][i]-best_WTG[t][i])*r*o, 0), P_WTG_MAX[t][i]);
                        }
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
            a_TGU = data.get("a_TGU");b_TGU = data.get("b_TGU");c_TGU = data.get("c_TGU");d_TGU = data.get("d_TGU");e_TGU = data.get("e_TGU");
            a_CHP = data.get("a_CHP");b_CHP = data.get("b_CHP");c_CHP = data.get("c_CHP");d_CHP = data.get("d_CHP");e_CHP = data.get("e_CHP");f_CHP = data.get("f_CHP");
            a_HOU = data.get("a_HOU");b_HOU = data.get("b_HOU");c_HOU = data.get("c_HOU");
            UR_TGU = data.get("UR_TGU");DR_TGU = data.get("DR_TGU");UR_CHP = data.get("UR_CHP");DR_CHP = data.get("DR_CHP");
            T_ref = data.get("T_ref");T_amb = data.get("T_amb");a_SPVP = data.get("a_SPVP");P_SPVP_Reted = data.get("P_SPVP_Rated");S_SPVP = data.get("S_SPVP");CC_SPVP = data.get("CC_SPVP");
            V_Wind_CUT_IN = data.get("V_Wind_CUT_IN");V_Wind_CUT_OFF = data.get("V_Wind_CUT_OFF");V_Wind_Rated = data.get("V_Wind_Rated");P_Wind_Rated = data.get("P_Wind_Rated");V_Wind_Real = data.get("V_Wind_Real");CC_WTG = data.get("CC_WTG");

            P_CHP_MAX = data.get("P_CHP_MAX");P_CHP_MIN = data.get("P_CHP_MIN");H_CHP_MIN = data.get("H_CHP_MIN");H_CHP_MAX = data.get("H_CHP_MAX");H_CHP_MID = data.get("H_CHP_MID");
            P_TGU_MAX = data.get("P_TGU_MAX");P_TGU_MIN = data.get("P_TGU_MIN");
            H_HOU_MIN = data.get("H_HOU_MIN");H_HOU_MAX = data.get("H_HOU_MAX");

            P_Need = data.get("P_Need");H_Need = data.get("H_Need");
            fitVal = new double[popSize];
            pbestFitness = 1e20;


            P_P2G_MAX = new double[]{30,20,30,20,5,20};
            P_CHP_CUR_MAX = new double[popSize][25][6];P_CHP_CUR_MIN = new double[popSize][25][6];
            P_SPVP_MAX = new double[25][3];P_WTG_MAX = new double[25][3];

            p_P2G = new double[popSize][25][6];v_P2G = new double[popSize][25][6];prep_P2G = new double[popSize][25][6];
            p_TGU = new double[popSize][25][13];v_TGU = new double[popSize][25][13];prep_TGU = new double[popSize][25][13];
            p_CHP = new double[popSize][25][6];v_p_CHP = new double[popSize][25][6];prep_CHP = new double[popSize][25][6];
            p_SPVP = new double[popSize][25][3];v_SPVP = new double[popSize][25][3];prep_SPVP = new double[popSize][25][3];
            p_WTG = new double[popSize][25][3];v_WTG = new double[popSize][25][3];prep_WTG = new double[popSize][25][3];
            h_CHP = new double[popSize][25][6];v_h_CHP = new double[popSize][25][6];preh_CHP = new double[popSize][25][6];
            h_HOU = new double[popSize][25][3];v_HOU = new double[popSize][25][3];preh_HOU = new double[popSize][25][3];

            best_TGU = new double[25][13];
            best_p_CHP = new double[25][6];
            best_h_CHP = new double[25][6];
            best_SPVP = new double[25][3];
            best_WTG = new double[25][3];
            best_HOU = new double[25][3];
            best_P2G = new double[25][6];
            windPower();solarPower();
            for(int s = 0; s < popSize; s++){
                System.out.println(s);
                boolean flag = true;
                do{
                    flag = true;
                    for(int t = 1; t <= 24 && flag; t++){
                        for(int i = 0; i < N_CHP; i++) h_CHP[s][t][i] = random.nextDouble(H_CHP_MAX[i]*0.95,H_CHP_MAX[i]);
                        for(int i = 0; i < N_HOU; i++) h_HOU[s][t][i] = random.nextDouble(H_HOU_MAX[i]*0.7,H_HOU_MAX[i]);
                        for(int i = 0; i < N_CHP; i++) p_P2G[s][t][i] = random.nextDouble(0,P_P2G_MAX[i]);
                        if(!constrainH(s,t)) flag = false;
                    }
                    if(!flag) continue;
                    CHP_FOR(s);
                    for(int t = 1; t <= 24 && flag; t++){
                        if(t >= 2){
                            for(int i = 0; i < N_TGU; i++) p_TGU[s][t][i] = random.nextDouble(Math.max(p_TGU[s][t-1][i] - DR_TGU[i],P_TGU_MIN[i]), Math.min(p_TGU[s][t-1][i] + UR_TGU[i],P_TGU_MAX[i]));
                            for(int i = 0; i < N_CHP; i++){
                                p_CHP[s][t][i] = random.nextDouble(Math.max(p_CHP[s][t-1][i] - DR_CHP[i],P_CHP_CUR_MIN[s][t][i]), Math.min(p_CHP[s][t-1][i] + UR_CHP[i],P_CHP_CUR_MAX[s][t][i]));
                            }
                        }else{
                            for(int i = 0; i < N_TGU; i++) p_TGU[s][t][i] = random.nextDouble(Math.max(P_TGU_MAX[i]*0.6, P_TGU_MIN[i]),P_TGU_MAX[i]);
                            for(int i = 0; i < N_CHP; i++) p_CHP[s][t][i] = random.nextDouble(P_CHP_CUR_MIN[s][t][i],P_CHP_CUR_MAX[s][t][i]);
                        }
                        for(int i = 0; i < N_SPVP; i++) p_SPVP[s][t][i] = P_SPVP_MAX[t][i] == 0 ? 0 : random.nextDouble(P_SPVP_MAX[t][i]*0.95,P_SPVP_MAX[t][i]);
                        for(int i = 0; i < N_WTG; i++) p_WTG[s][t][i] = P_WTG_MAX[t][i] == 0 ? 0 : random.nextDouble(P_WTG_MAX[t][i]*0.95,P_WTG_MAX[t][i]);
                        if(!constrainP(s,t)) flag = false;
                    }
                }while(!flag);

                for(int t = 1; t <= 24; t++){
                    for(int i = 0; i < N_TGU; i++){
                        prep_TGU[s][t][i] = p_TGU[s][t][i];
                        v_TGU[s][t][i] = random.nextDouble(P_TGU_MIN[i]/10, P_TGU_MIN[i]/10 + (P_TGU_MAX[i]-P_TGU_MIN[i])/10);
                    }
                    for(int i = 0; i < N_CHP; i++){
                        prep_CHP[s][t][i] = p_CHP[s][t][i];
                        v_p_CHP[s][t][i] = random.nextDouble(P_CHP_CUR_MIN[s][t][i]/10, P_CHP_CUR_MIN[s][t][i]/10 + (P_CHP_CUR_MAX[s][t][i]-P_CHP_CUR_MIN[s][t][i])/10);
                        preh_CHP[s][t][i] = h_CHP[s][t][i];
                        v_h_CHP[s][t][i] = random.nextDouble(H_CHP_MIN[i]/10, H_CHP_MIN[i]/10 + (H_CHP_MAX[i]-H_CHP_MIN[i])/10);
                        prep_P2G[s][t][i] = p_P2G[s][t][i];
                        v_P2G[s][t][i] = random.nextDouble(0, P_P2G_MAX[i]/10);
                    }
                    for(int i = 0; i < N_HOU; i++){
                        preh_HOU[s][t][i] = h_HOU[s][t][i];
                        v_HOU[s][t][i] = random.nextDouble(H_HOU_MIN[i]/10, H_HOU_MIN[i]/10 + (H_HOU_MAX[i]-H_HOU_MIN[i])/10);
                    }
                    for(int i = 0; i < N_SPVP; i++){
                        prep_SPVP[s][t][i] = p_SPVP[s][t][i];
                        v_SPVP[s][t][i] = random.nextDouble(0, P_SPVP_MAX[t][i]/10+0.000001);
                    }
                    for(int i = 0; i < N_WTG; i++){
                        prep_WTG[s][t][i] = p_WTG[s][t][i];
                        v_WTG[s][t][i] = random.nextDouble(0, P_WTG_MAX[t][i]/10+0.000001);
                    }
                }
            }
            fitness();
        }
        public static void myQuicksort(int low, int high){
            if(low < high){
                int l = low, r = high;
                double tmp = fitVal[l];
                double[][] tmp1 = p_TGU[l], tmp2 = p_CHP[l], tmp3 = h_CHP[l], tmp4 = h_HOU[l], tmp5 = p_SPVP[l],
                        tmp6 = p_WTG[l], tmp7 = p_P2G[l], tmp8 = v_TGU[l], tmp9 = v_p_CHP[l], tmp10 = v_h_CHP[l],
                        tmp11 = v_HOU[l], tmp12 = v_SPVP[l], tmp13 = v_WTG[l], tmp14 = v_P2G[l];
                while(l < r){
                    while(l < r && fitVal[r] > tmp) r--;
                    if(l < r){
                        fitVal[l] = fitVal[r];
                        p_TGU[l] = p_TGU[r]; p_CHP[l] = p_CHP[r]; h_CHP[l] = h_CHP[r]; h_HOU[l] = h_HOU[r]; p_SPVP[l] = p_SPVP[r]; p_WTG[l] = p_WTG[r];
                        p_P2G[l] = p_P2G[r]; v_TGU[l] = v_TGU[r]; v_p_CHP[l] = v_p_CHP[r]; v_h_CHP[l] = v_h_CHP[r]; v_HOU[l] = v_HOU[r];
                        v_SPVP[l] = v_SPVP[r]; v_WTG[l] = v_WTG[r]; v_P2G[l] = v_P2G[r];
                        l++;
                    }
                    while(l < r && fitVal[l] < tmp) l++;
                    if(l < r){
                        fitVal[r] = fitVal[l];
                        p_TGU[r] = p_TGU[l]; p_CHP[r] = p_CHP[l]; h_CHP[r] = h_CHP[l]; h_HOU[r] = h_HOU[l]; p_SPVP[r] = p_SPVP[l]; p_WTG[r] = p_WTG[l];
                        p_P2G[r] = p_P2G[l]; v_TGU[r] = v_TGU[l]; v_p_CHP[r] = v_p_CHP[l]; v_h_CHP[r] = v_h_CHP[l]; v_HOU[r] = v_HOU[l];
                        v_SPVP[r] = v_SPVP[l]; v_WTG[r] = v_WTG[l]; v_P2G[r] = v_P2G[l];
                        r--;
                    }
                }
                fitVal[l] = tmp;
                p_TGU[l] = tmp1; p_CHP[l] = tmp2; h_CHP[l] = tmp3; h_HOU[l] = tmp4; p_SPVP[l] = tmp5;
                p_WTG[l] = tmp6; p_P2G[l] = tmp7; v_TGU[l] = tmp8; v_p_CHP[l] = tmp9; v_h_CHP[l] = tmp10;
                v_HOU[l] = tmp11; v_SPVP[l] = tmp12; v_WTG[l] = tmp13; v_P2G[l] = tmp14;
                myQuicksort(low, l-1);
                myQuicksort(l+1, high);
            }
        }
        public void start(String filename) {
            FileOutputStream foss = null;
            try {
                foss = new FileOutputStream("file2.txt",true);
            }catch (Exception ignored){}
            init();
            for(int its = 1; its <= T; its++){
                double xita = random.nextDouble(0,1)*2*Math.PI;
                myQuicksort(0,popSize-1);
                updateV(its,xita);
                if(xita <= Math.PI) updateStageOne(4*random.nextDouble(0,1)-2, 3*random.nextDouble(0,1)-1.5, 2*random.nextDouble(0,1)-1);
                else updateStageTwo(its);
                fitness();
                System.out.println("iterations : " + its + " / " + T + "  | minFit : " + pbestFitness);
                try {
                    assert foss != null;
                    foss.write((String.valueOf(pbestFitness) + '\n').getBytes());
                }catch (Exception ignored){}
            }
            double[] pb = new double[25], hb = new double[25];
            for(int t = 1; t <= 24; t++){
                for(int i = 0; i < 13; i++) pb[t] += best_TGU[t][i];
                for(int i = 0; i < 6; i++) pb[t] += (best_p_CHP[t][i] + best_P2G[t][i]*0.4*0.55);
                for(int i = 0; i < 3; i++) pb[t] += best_SPVP[t][i];
                for(int i = 0; i < 3; i++) pb[t] += best_WTG[t][i];
                for(int i = 0; i < 6; i++) hb[t] += best_h_CHP[t][i];
                for(int i = 0; i < 3; i++) hb[t] += best_HOU[t][i];
                pb[t] -= P_Need[t]*1.1;
                hb[t] -= H_Need[t]*1.15;
            }


            try {
                assert foss != null;
                foss.close();
                FileOutputStream fos = new FileOutputStream(filename,true);
                String s = String.valueOf(pbestFitness) + " |  " + Arrays.deepToString(best_TGU) + "   " + Arrays.deepToString(best_p_CHP)
                        + "   "+ Arrays.deepToString(best_SPVP) + "   "+ Arrays.deepToString(best_WTG)
                        + "   " + Arrays.deepToString(best_P2G) + "   "+ Arrays.deepToString(best_h_CHP)
                        + "   " + Arrays.deepToString(best_HOU) + '\n';
                fos.write(s.getBytes());
                s = Arrays.toString(pb) + "      " + Arrays.toString(hb) + '\n' + '\n';
                fos.write(s.getBytes());
                fos.close();
            }catch (Exception ignored){}
        }
    }

    public static void main(String[] args) {
        for(int eachtime = 0 ; eachtime < 30; eachtime++){
            SGA sga = new SGA(60, 5000, 0.8, 20, 0);
            String filename = "file1.txt";
            sga.start(filename);
        }
    }
}


