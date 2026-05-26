`timescale 1ns / 1ps

module N_FIR
#(
    parameter WI = 2, WF = 10,
    parameter WIC = 1, WFC = 15,
    parameter N_TAP = 6,
    parameter WIO = N_TAP-1 + (WI+WIC), WFO = WF+WFC
)(
    input CLK, RST,
    input [(N_TAP * (WIC+WFC))-1:0] flat_h,
    input signed [(WI+WF)-1:0] xn,
    output signed [(WIO+WFO)-1:0] yn
);
    reg [(WIC+WFC)-1:0] h_mem [0:N_TAP-1];
    integer i;
    
    always @* begin
        for(i = 0; i < N_TAP; i = i + 1) begin
            h_mem[i] = flat_h[(i*(WIC+WFC)) +: WIC+WFC];
        end
    end
    
    reg signed [(WI+WF)-1:0] d_ff [0:N_TAP-1];                  // 2D Delayed register variable
    wire signed [(WI+WF+WIC+WFC)-1:0] mult_result [0:N_TAP-1];  // 2D Multiplier result variable
    wire signed [(N_TAP-1) + (WI+WF+WIC+WFC+1)-1:0] add_result [0:N_TAP-2];   // 2D Adder result variable

    
    always @(posedge CLK) begin
        if (RST) begin                                      // If RST is active high:
            for(i = 0; i < N_TAP; i = i + 1) begin
                d_ff[i] <= 0;                               // Set the delayed registers to zero
            end
        end    
        else begin                                          // If RST is active low: 
            d_ff[0] <= xn;                                  // Assign first element of delayed register to input
            for(i = 1; i < N_TAP; i = i + 1) begin
                d_ff[i] <= d_ff[i-1];                   // Shift each delay register elment to the right
            end
        end
    end
    
    generate
        genvar k;
        for(k = 0; k < N_TAP; k = k + 1) begin // MAC Units
            // Multiplier
            FP_multiplier #(.WI1(WI),.WF1(WF),.WI2(WIC),.WF2(WFC),.WIO(WI+WIC),.WFO(WF+WFC))
                        FM (.in1(d_ff[k]),
                            .in2(h_mem[k]),
                            .out(mult_result[k]));  
            
            // Adder                
            if (k == 0) begin
            FP_adder #(.WI1(WI+WIC),.WF1(WF+WFC),.WI2(WI+WIC),.WF2(WF+WFC),.WIO(WI+WIC+1),.WFO(WF+WFC))
                        FA (.in1(0),
                            .in2(mult_result[k]),
                            .out(add_result[k]));  
            end
            else begin
            FP_adder #(.WI1(WI+WIC+(k-1)),.WF1(WF+WFC),.WI2(WI+WIC),.WF2(WF+WFC),.WIO(WI+WIC+(k-1)+1),.WFO(WF+WFC))
                        FA (.in1(add_result[k-1]),
                            .in2(mult_result[k]),
                            .out(add_result[k])); 
            end 
        end          
    endgenerate
    
    assign yn = add_result[N_TAP-2];
       
endmodule
