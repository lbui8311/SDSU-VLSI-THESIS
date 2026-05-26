`timescale 1ns / 1ps

module FP_Mitchell_Approx
#(
    parameter WI = 5, WF = 11,
    parameter WIO = 2, WFO = 14
)(
    input RST,
    input [(WI+WF)-1:0] x,                // Q(5.11)
    output reg signed [(WIO+WFO)-1:0] y       // Q(2.14)
);   

    reg [$clog2(WI+WF)-1:0] k;
    reg [(WI+WF)-1:0] temp_val;
    reg valid;
    
    // Leading one detector
     priorty_encoder #(.N(WI+WF))
             ENC0 (.in(x),.code(k),.valid(valid)); 

    always @* begin
        // Mantissa Extration => Shift leading 1 to MSB
        //temp_val = x << ((WI+WF-1)-k);
        temp_val = {{WF{1'b0}}, x};
        temp_val = (temp_val << WF) >> k;
        
        if (x == 0) y = 0; // log2(0) is mathematically undefined
        else y = {k, temp_val[WF-1:0]};
    end
    
endmodule
