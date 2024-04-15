package EconomicDispatch;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class GetMean {
    public static void main(String[] args) throws IOException {
        String filename = "PSOData1.txt";
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String s = br.readLine();
        double[] nums = new double[30];
        int numsSize = 0;
        double mean = 0;
        double std = 0;
        while(s != null){
            nums[numsSize++] = Double.parseDouble(s.substring(0,s.indexOf("|")-1));
            mean += nums[numsSize-1];
            s = br.readLine();
        }
        mean = mean / 30;

        for(int i = 0; i < 30; i++){
            std += Math.pow(mean - nums[i],2);
        }
        std = Math.sqrt(std);
        System.out.println("mean = " + mean + "   std = " + std);
    }
}
