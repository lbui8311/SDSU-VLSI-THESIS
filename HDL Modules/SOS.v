`timescale 1ns / 1ps

module SOS // DF2 Transposed
#(
    parameter WI = 2, WF = 10,
    parameter WIC = 2, WFC = 14,
    //parameter WI_MAC = WI + WIC + 3, WF_MAC = WF + WFC,
    parameter WI_MAC = WI + WIC + 3, WF_MAC = WF+(WFC * 2 + 2),
    parameter WIO = WI_MAC, WFO = WF_MAC
)(
    input CLK, RST,
    input signed [(WI+WF)-1:0] xn,
    input signed [(WIC+WFC)-1:0] a1, a2, b1, b2,
    output signed [(WIO+WFO)-1:0] yn
);
    
    // Delay Line
    reg signed [(WI_MAC+WF_MAC)-1:0] v_n1, v_n2;
    
    // Computation outputs
    wire signed [(WI_MAC+WF_MAC)-1:0] a1_product, a2_product, b1_product, b2_product;
    wire signed [(WI_MAC+WF_MAC)-1:0] a_sum, b_sum, v_n;
       
    // Delay Sequential Logic
    always @(posedge CLK) begin
        if (RST) begin
            v_n1 <= 0; v_n2 <= 0;
        end
        else begin
            v_n1 <= v_n; v_n2 <= v_n1;
        end
    end
    
    // Feedback Path (left)
    // v(n) = x(n) + (-a1*v(n-1)) + (-a2*v(n-2))
    FP_multiplier #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WIC),.WF2(WFC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FM_a1 (.in1(v_n1),
                               .in2(a1),
                               .out(a1_product)); 
    FP_multiplier #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WIC),.WF2(WFC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FM_a2 (.in1(v_n2),
                               .in2(a2),
                               .out(a2_product));   
                                                  
    FP_adder #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WI_MAC),.WF2(WF_MAC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FA_sa (.in1(a1_product),
                               .in2(a2_product),
                               .out(a_sum)); 
                            
    FP_adder #(.WI1(WI),.WF1(WF),.WI2(WI_MAC),.WF2(WF_MAC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FA_vn (.in1(xn),
                               .in2(a_sum),
                               .out(v_n));                         
    
    // Feedforward Path (right)
    // y(n) = b1*v(n-1) + b2*v(n-2)
    FP_multiplier #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WIC),.WF2(WFC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FM_b1 (.in1(v_n1),
                               .in2(b1),
                               .out(b1_product)); 
    FP_multiplier #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WIC),.WF2(WFC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FM_b2 (.in1(v_n2),
                               .in2(b2),
                               .out(b2_product));   
                                                  
    FP_adder #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WI_MAC),.WF2(WF_MAC),.WIO(WI_MAC),.WFO(WF_MAC))
                        FA_sb (.in1(b1_product),
                               .in2(b2_product),
                               .out(b_sum)); 
                            
    FP_adder #(.WI1(WI_MAC),.WF1(WF_MAC),.WI2(WI_MAC),.WF2(WF_MAC),.WIO(WIO),.WFO(WFO))
                        FA_yn (.in1(v_n),
                               .in2(b_sum),
                               .out(yn)); 
                               
endmodule