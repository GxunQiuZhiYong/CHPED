package EconomicDispatch;
public class GetIntegrals {
    public static double get(String mode,double a, double b) {
        if(a == b) return 0;
        double accuracy;   //a为区间下限，b为区间上限，accuracy为精度
        accuracy = 1e-3;//精度
        double h;   //步长
        h = b - a;
        double T1, T2;   //T1：二分前的梯形法积分值；T2：二分后的梯形法积分值；
        T1 = h / 2 * (1 + f(b, mode,a,b));
        T2 = 0;
        int k = 1;  //记录二分次数
        double S1 = 0, S2;   //对T1与T2加权平均，得辛普森积分值
        double C1=0, C2;   //对S1与S2加权平均，得柯特斯积分值
        double R1 = 0, R2;   //对C1，C2加权平均，得龙贝格积分值
        while (true){
            double sum = 0;  //各分点的函数值和
            double x = a + h / 2;   //分点值
            while (x < b)   //在区间上限范围内求各分点的函数值和
            {
                sum += f(x, mode,a,b);
                x += h;
            }
            T2 = T1 / 2 + h / 2 * sum;    //计算梯形序列得下一个二分结果
            S2 = T2 + 1.0 / 3 * (T2 - T1);   //线性组合外推值simpson
            if (k == 1)   //至少外推2次得出S1，S2
            {
                k++;
                h /= 2;
                T1 = T2;
                S1 = S2;
            }
            else
            {
                C2 = S2 + 1.0 / 15 * (S2 - S1);  //线性组合外推值Cotes
                if (k == 2)   //至少外推3次得出C1，C2
                {
                    C1 = C2;
                    k++;
                    h /= 2;
                    T1 = T2;
                    S1 = S2;
                }
                else
                {
                    R2 = C2 + 1.0 / 63 * (C2 - C1);   //线性组合外推至Romberg

                    if (k == 3)   //至少外推4次得出R1,R2
                    {
                        R1 = R2;
                        C1 = C2;
                        k++;
                        h /= 2;
                        T1 = T2;
                        S1 = S2;
                    }
                    else if (Math.abs(R2 - R1) >= accuracy)   //精度仍然不符合要求，继续二分步长、继续外推
                    {
                        R1 = R2;
                        C1 = C2;
                        k = k + 1;
                        h = h / 2;
                        T1 = T2;
                        S1 = S2;
                    }
                    else   //精度符合要求，修改flag为0，跳出while循环
                    {
                        return R2;
                    }
                }
            }
        }
    }
    static double f(double x, String mode, double a, double b)   //自定义被积函数
    {
        double result;
        if(mode.equals("SPVP_Oe")) result = (b-x)*2.95175/x*Math.pow((x/11.7714),2.95175)*Math.pow(Math.E, -Math.pow(x/11.7714,2.95175));
        else if (mode.equals("SPVP_Ue")) result = (x-a)*2.95175/x*Math.pow((x/11.7714),2.95175)*Math.pow(Math.E, -Math.pow(x/11.7714,2.95175));
        else if (mode.equals("WTG_Oe")) result = (b-x)*0.522/(Math.pow(Math.PI*2, 0.5)*x)*Math.pow(Math.E,-Math.pow(Math.log(x/5.278)/Math.log(Math.E),2)/(0.5*0.522*0.522));
        else result = (x-a)*0.522/(Math.pow(Math.PI*2, 0.5)*x)*Math.pow(Math.E,-Math.pow(Math.log(x/5.278)/Math.log(Math.E),2)/(0.5*0.522*0.522));
        return result;
    }

    public static void main(String[] args) {
        System.out.println(get("SPVP_Oe",10,80));
        System.out.println(get("SPVP_Ue",0,175));
        System.out.println(get("WTG_Oe",0,175));
        System.out.println(get("WTG_Ue",0,175));
    }
}
