`timescale 1ns / 1ps

module FP_Z_Score
#(
    parameter k = 2,
    parameter WI = 5, WF = 11, WL = WI + WF,
    parameter WIO = 2, WFO = 14, 
    parameter WL_coeff = 32, WL_LUT = 6
)(
    input RST_n, CLK,
    input signed [(WI+WF)-1:0] x_in,                // Q(5.11)
    input signed [(WI+WF)-1:0] alpha, epsilon, 
    output signed [(WIO+WFO)-1:0] z_out,       // Q(2.14)
    output z_valid,
    output overflow
);     

    // Input pre-scale
    wire signed [WL-1:0] x_scaled = x_in >>> k;

    // State registers
    reg signed [WL-1:0] mean_reg;    // mean out
    reg signed [WL-1:0] var_reg;     // var out
    reg signed [WL-1:0] alpha_x_reg; // x_scaled * a 
    reg signed [WL-1:0] x_sqr_reg;   // (x_scaled - mean out)^2
    reg signed [WL-1:0] one_a_reg;   // 1-a
    reg        [4:0] w_count;        // Window count
    
    // 1 - alpha
    wire signed [WL-1:0] one_a_w;
    
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
             FA0 (.in1(1),.in2(-alpha),.out(one_a_w));  
             
    // Mean tracking = a*x_t + (1 - a)*u_t-1
    wire signed [WL-1:0] a_x_t, one_a_u_t, a_x_one_a_u_t;
    
    FP_mult #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
                  FP0 (.in1(alpha),.in2(x_scaled),.out(a_x_t)); 
                  
    FP_mult #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
                  FP1 (.in1(one_a_w),.in2(mean_reg),.out(one_a_u_t)); 
    
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
             FA1 (.in1(a_x_t),.in2(one_a_u_t),.out(a_x_one_a_u_t)); 
             
    // Deviation = x_s - mean_reg
    wire signed [WL-1:0] dev_w;
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
             FA2 (.in1(x_scaled),.in2(-mean_reg),.out(dev_w)); 
    
    // Variance tracking = a*(x_t - u_t)^2 + (1-a)*v_t-1        
    wire signed [WL-1:0] dev_sqr, dev_sqr_a, var_tmp, var_nxt;           
    
    FP_mult #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
                  FP2 (.in1(dev_w),.in2(dev_w),.out(dev_sqr)); 
                  
    FP_mult #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
                  FP3 (.in1(alpha),.in2(dev_sqr),.out(dev_sqr_a)); 
    
    FP_mult #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
                  FP4 (.in1(one_a_w),.in2(var_reg),.out(var_tmp)); 
                            
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
             FA3 (.in1(dev_sqr_a),.in2(var_tmp),.out(var_nxt));  
     
    // Adding epsilon to v(t)
    wire signed [WL-1:0] var_eps;
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
             FA4 (.in1(var_nxt),.in2(epsilon),.out(var_eps));  
             
    // Log2(x) Variance
    wire signed [WL-1:0] var_log;
    FP_Mitchell_Approx #(.WI(WI),.WF(WF),.WIO(WI),.WFO(WF))
             LG0 (.x(var_eps),.y(var_log));           
    
    // >>> 1 for sqrt
    wire signed [WL-1:0] lg_shft = var_log >>> 1;
    
    // |x|
    wire signed [WL-1:0] abs_dev = dev_w[WL-1] ? - dev_w : dev_w;
    
    // Log(2) of |x|
    wire signed [WL-1:0] dev_log;
    FP_Mitchell_Approx #(.WI(WI),.WF(WF),.WIO(WI),.WFO(WF))
             LG1 (.x(dev_w),.y(dev_log));  
    
    // var_log >>> 1 + dev_log of |x|
    wire signed [WL-1:0] diff_var_dev;
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI),.WFO(WF))
             FA5 (.in1(lg_shft),.in2(dev_log),.out(diff_var_dev));  
   
    // Inverse log
    wire signed [WL-1:0] z_inv;
    
    FP_Inv_Mitchell #(.WI(WI),.WF(WF),.WIO(WI),.WFO(WF))
             LG2 (.x(diff_var_dev),.y(z_inv));    
             
    wire signed [WL-1:0] z_temp = z_inv >>> k;
    
    always @(posedge CLK or negedge RST_n) begin
        if (!RST_n) begin
            mean_reg <= {WL{1'b0}};
            var_reg <= {WL{1'b0}};
            alpha_x_reg <= {WL{1'b0}};
            x_sqr_reg <= {WL{1'b0}};
            one_a_reg <= {WL{1'b0}};
        end else begin
            alpha_x_reg <= a_x_t;
            one_a_reg <= one_a_w;
            x_sqr_reg <= dev_w;
            
            if (w_count == 5'd0) begin // if first value
                mean_reg <= x_in;
                var_reg <= {WL{1'b0}};
            end else begin
                mean_reg <= a_x_one_a_u_t;
                var_reg <= var_nxt;
            end
            
            if (w_count != {5{1'b1}}) w_count <= w_count + 5'd1;  
        end
    end                         
endmodule

