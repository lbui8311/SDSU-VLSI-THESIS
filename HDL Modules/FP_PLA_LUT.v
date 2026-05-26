`timescale 1ns / 1ps

module FP_PLA_LUT 
#(
    parameter WI = 5, WF = 11,
    parameter WIO = 2, WFO = 14, 
    parameter WL_coeff = 32, WL_LUT = 6
)(
    input RST,
    input [(WI+WF)-1:0] x,                // Q(5.11)
    output signed [(WIO+WFO)-1:0] y       // Q(2.14)
);     
    
    reg signed [2*(WI+WF)-1:0] rom_LUT [0:(2**WL_LUT)-1]; // 2-D Lookup Table
    reg signed [WI+WF-1:0] a, b;                      // Coefficient values, each 16-bit
    reg [$clog2(WL_LUT)-1:0] segment_index;        // Output of Segment Decoder 6-bit
    wire [2*(WI+WF)-1:0] temp_index;                  // placeholder for index calculations
    reg [2*(WI+WF)-1:0] select_row;                   // row values of coefficients 'a' & 'b'
    wire signed [2*(WI+WF)-1:0] mult_result;
    integer i;
    reg signed [WL_coeff-1:0] a_coeff [0:(2**WL_LUT)-1];
    reg signed [WL_coeff-1:0] b_coeff [0:(2**WL_LUT)-1];

    // Thresholds
    localparam T1 = 16'b0000110010001111; // 1.57
    localparam T2 = 16'b0001100100011111; // 3.14
    localparam T3 = 16'b0010010110101110; // 4.71
    
    initial begin
        $readmemh("C:\Mac\Home\Downloads\Thesis\Thesis.srcs\sources_1\LUT_a.mem", a_coeff);
        $readmemh("C:\Mac\Home\Downloads\Thesis\Thesis.srcs\sources_1\LUT_b.mem", b_coeff);
    end
    
    initial begin // LUT Q(5.11)                             a        &      b           from MATLAB
            rom_LUT[0] = 32'b11111010101111000000100100111001; // -0.658203125000000 & -0.658203125000000
            rom_LUT[1] = 32'b11111010101011000000011101111000; // 1.152832031250000 &  -0.666015625000000
            rom_LUT[2] = 32'b00000101010101001110010111111011; // 0.933593750000000 & 0.666015625000000
            rom_LUT[3] = 32'b00000101010001001110100000100110; // -3.252441406250000 & 0.658203125000000
            for(i = 4; i < 64; i = i + 1) rom_LUT[i] = 32'd0;
    end

        
    always @* begin
        if (RST) begin // If RST is active high, clear index, row values and coefficient placeholders.
            segment_index = 0; 
            select_row = 0;
            a = 0; b = 0; 
        end
        else begin // else RST is active low
            //segment_index = temp_index[2*(WI+WF)-1:2*WF]; // concatenate the integer part from FP0 to obtain index Q(10.22)
            if (x < T1) segment_index = 3'd0;
            else if (x < T2) segment_index = 3'd1;
            else if (x < T3) segment_index = 3'd2;
            else segment_index = 3'd3;
            select_row = rom_LUT[segment_index];          // map select_row with corresponding index
            a = select_row[2*(WI+WF)-1:WI+WF];            // concatenate first half for a
            b = select_row[(WI+WF)-1:0];                  // concatenate second half for b
        end 
    end 
    
    // mult_result = ax
    FP_mult #(.WI1(WI),.WF1(WF),.WI2(WI),.WF2(WF),.WIO(WI+WI),.WFO(WF+WF))
                  FP1 (.in1(x),.in2(a),.out(mult_result)); 
    
    // y = mult_result + b
    FP_adder #(.WI1(WI+WI),.WF1(WF+WF),.WI2(WI),.WF2(WF),.WIO(WIO),.WFO(WFO))
             FA0 (.in1(mult_result),.in2(b),.out(y));             
endmodule
