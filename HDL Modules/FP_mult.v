`timescale 1ns / 1ps

`define WI_max (WI1 + WI2)
`define WF_max (WF1 + WF2)

module FP_mult #(parameter WI1 = 5, WF1 = 4,            // (5.4)
                       WI2 = 7, WF2 = 3,                // (7.3) 
                       WIO = 12, WFO = 7)          
                      (input signed [WI1+WF1-1:0] in1,
                       input signed [WI2+WF2-1:0] in2,
                       output signed [WIO+WFO-1:0] out); // (12.7)
    
    // Overflow flag
    wire OVF;
    
    //Temporary place holders for sign extend and zero pad
    reg signed [(WI1+WF1+WI2+WF2)-1:0] temp_out;
    reg [(WI1+WI2)-1:0] temp_int;
    reg [(WF1+WF2)-1:0] temp_frc;
    reg [WIO-1:0] int_out;
    reg [WFO-1:0] frc_out;
    
    // Find the current max output possible
    localparam signed [WIO-1:0] max_output = {1'b0, {(WIO-1){1'b1}}};
    
    always @* begin  
          temp_out = in1 * in2;
          
          temp_int = temp_out[`WI_max + `WF_max-1-:`WI_max]; // extract the integer part (left of radix)
          temp_frc = temp_out[`WF_max:0];  // extract the fractional part (right of radix)
                    
          if(WIO < `WI_max) int_out = {temp_int[`WI_max-1], temp_int[`WI_max-2:0]}; // truncate
          if(WIO >= `WI_max) int_out = {{(WIO-`WI_max){temp_int[`WI_max-1]}}, temp_int}; // sign extend  
          
          if(WFO < `WF_max) frc_out = temp_frc[`WF_max-1-:WFO]; // truncate
          if(WFO >= `WF_max) frc_out = {temp_frc, {(WFO-`WF_max){1'b0}}}; // zero append                   
            
    end 
       assign OVF = (temp_int > max_output) ? 1'b1 : 1'b0; 
       assign out = {int_out, frc_out};                                      
endmodule