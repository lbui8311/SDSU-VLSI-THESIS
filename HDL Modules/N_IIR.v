`timescale 1ns / 1ps

module N_IIR
#(
    parameter WI = 2, WF = 10,
    parameter WIC = 2, WFC = 14,
    parameter N = 2,
    parameter WIO = WI + WIC + 3, WFO = WF+(WFC * 2 + 2)
    //parameter WIO = WI + WIC + 3, WFO = WF + WFC
)(
    input CLK, RST,
    input [(N * 2 * (WIC+WFC))-1:0] flat_g,
    input [(N * (WIC+WFC))-1:0] a1,a2,b1,b2,
    input signed [(WI+WF)-1:0] xn,
    output signed [(WIO+WFO)-1:0] yn
);
    reg signed [2 * (WIC+WFC)-1:0] g_mem [0:N-1];
    reg signed [(WIC+WFC)-1:0] a1_mem [0:N-1];
    reg signed [(WIC+WFC)-1:0] a2_mem [0:N-1];
    reg signed [(WIC+WFC)-1:0] b1_mem [0:N-1];
    reg signed [(WIC+WFC)-1:0] b2_mem [0:N-1];
    integer i, j;
    
    // Unrolling the gains
    always @* begin
        for(i = 0; i < N; i = i + 1) begin
            g_mem[i] = flat_g[(i*(2 * (WIC+WFC))) +: 2 * (WIC+WFC)];
        end
    end
    
    // Unrolling the coefficients
    always @* begin
        for(j = 0; j < N; j = j + 1) begin
            a1_mem[j] = a1[(j*(WIC+WFC)) +: WIC+WFC];
            a2_mem[j] = a2[(j*(WIC+WFC)) +: WIC+WFC];
            b1_mem[j] = b1[(j*(WIC+WFC)) +: WIC+WFC];
            b2_mem[j] = b2[(j*(WIC+WFC)) +: WIC+WFC];
        end
    end
    
    wire signed [2 * (WI+WF)-1:0] g_products [0:N-1];
    wire signed [(WIO+WFO)-1:0] SOS_out [0:N-1];
    wire signed [(WIO+WFO)-1:0] stage_input [0:N-1];
    
    generate
        genvar k;
        for(k = 0; k < N; k = k + 1) begin // Scalling down the input to each SOS
            // if k = 0, input is xn. If k > 0, input is output of previous SOS
            assign stage_input[k] = (k == 0) ? xn : SOS_out[k-1];
            FP_multiplier #(.WI1(WI),.WF1(WF),.WI2(WIC),.WF2(WFC * 2 + 2),.WIO(WI+WIC),.WFO(WF+(WFC * 2 + 2)))
                        FM (.in1(stage_input[k]),
                            .in2(g_mem[k]),
                            .out(g_products[k]));  
//            fp_multiplier #(.WI1(WI),.WF1(WF),.WI2(WIC),.WF2(WFC * 2 + 2),.WIO(WI+WIC),.WFO(WF+(WFC * 2 + 2)))
//                        FM (.in1(stage_input[k]),
//                            .in2(g_mem[k]),
//                            .out(full_g_products[k]));
                     
//            assign g_products[k] = full_g_products[k][DROP_WL +: (WI+WF)];
                            
            SOS #(.WI(WI),.WF(WF),.WIC(WIC),.WFC(WFC),.WI_MAC(WI+WIC + 3),.WF_MAC(WF+(WFC * 2 + 2)),.WIO(WI),.WFO(WF))
                        SOS_N (.CLK(CLK),.RST(RST),
                               .xn(g_products[k]),
                               .a1(a1_mem[k]),.a2(a2_mem[k]),.b1(b1_mem[k]),.b2(b2_mem[k]),
                               .yn(SOS_out[k])); 
        end
    endgenerate 
   
    assign yn = SOS_out[N-1];    
     
endmodule
