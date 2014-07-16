//Angela Fan

//global variable for data
var data;

var file;

//margins and bounding boxes for the graph
var bb_linegraph;

bb_linegraph = {
	h: 300,
	w: 720,
	margin: {
		top: 50,
		right: 50,
		bottom: 20,
		left: 50
	}
};

svg_linegraph = d3.select("#linegraph").append("svg")
	.attr("class", "linegraph_svg")
	.attr({
		width: bb_linegraph.w + bb_linegraph.margin.left + bb_linegraph.margin.right + 400,
		height: bb_linegraph.h + bb_linegraph.margin.bottom + bb_linegraph.margin.top
	});

var linegraph = svg_linegraph.append("g")
	.attr("class", "linegraph")
	.attr("transform", "translate(" + (bb_linegraph.margin.left-40) + "," + bb_linegraph.margin.top + ")");

//tip call
var graph_tip = d3.tip()
	.attr("class", "d3-tip")
	.offset([0,0]);

svg_linegraph.call(graph_tip);

//function calls
draw_linegraph();


//linegraph globals
var linegraph_x, linegraph_y, linegraph_color, linegraph_xAxis, linegraph_yAxis;
var linegraph_xdomain, linegraph_ydomain;

function draw_linegraph() {

	linegraph_x = d3.scale.linear()
		.range([bb_linegraph.w, 0]);

	linegraph_y = d3.scale.linear()
		.range([bb_linegraph.h, 0]);

	linegraph_color = d3.scale.ordinal()
		.domain([0,1])
		.range(["#868686", "#b30000"]);

	linegraph_xAxis = d3.svg.axis()
		.scale(linegraph_x)
		.orient("bottom")
		//get rid of weird default end tick
		.outerTickSize([0]);

	linegraph_yAxis = d3.svg.axis()
		.scale(linegraph_y)
		.orient("right")
		//get rid of weird default end tick
		.outerTickSize([0]);

	//add the clip path for the brush
	linegraph.append("defs")
		.append("clipPath")
		.attr("transform", "translate(-4.5,-4.5)")
		.attr("id", "linegraph_clip")
		.append("rect")
		.attr("width", bb_linegraph.w+4.5)
		.attr("height", bb_linegraph.h + 4.5);

	linegraph.append("g")
		.attr("class", "linegraph_brush");

	linegraph.append("g")
		.attr("class", "x axis")
		.attr("transform", "translate(0," + bb_linegraph.h + ")")
		.call(linegraph_xAxis)
		.append("text")
		.attr("class", "label")
		.attr("y", 38)
		.attr("x", bb_linegraph.w/2 + 50)
		.attr("dy", "-.71em")
		.style("text-anchor", "end")
		.text("Enhancers Ranked by Signal");

	linegraph.append("g")
		.attr("class", "y axis")
		.attr("transform", "translate(" + bb_linegraph.w + ",0)")
		.call(linegraph_yAxis)
		.append("text")
		.attr("class", "label")
		.attr("transform", "rotate(-90)")
		.attr("y", 6)
		.attr("dy", "6.5em")
		.attr("x", -115)
		.style("text-anchor", "end")
		.text("Enhancer Signal");

	//titles 
	linegraph.append("text")
		.attr("y", -38)
		.attr("font-size", "16px")
		.attr("font-weight", "bold")
		.attr("text-anchor", "middle")
		.attr("x", bb_linegraph.w/2 + 30)
		.text("Hockey Stick");

	linegraph.append("text")
		.attr("y", -21)
		.attr("text-anchor", "middle")
		.attr("x", bb_linegraph.w/2 + 30)
		.style("font-size", "12px")
		.text("Select a region to zoom in. Click 'Clear Zoom' to zoom out.");

	linegraph.append("text")
		.attr("y", -8)
		.attr("text-anchor", "middle")
		.attr("x", bb_linegraph.w/2 + 30)
		.style("font-size", "12px")
		.text("Mouseover a circle more information, click to bring up a PDF, double-click to add to the table.");

	$(".dropdown-menu li a").click(function () {
		file = $(this).text()
		update_linegraph(file);

		d3.select(".linegraph_clear_button").remove();

		output_name = file.split("_HOCKEY");

		d3.select(".page_title").remove();

		d3.select("#page_title")
			.append("text")
			.attr("y", 0)
			.attr("text-anchor", "middle")
			.attr("x", 500)
			.text(file.split("_HOCKEY")[0])
			.attr("class", "page_title");
	});

	d3.csv("/Documents/Bradner_work/hockey-sticks/lookup_table.csv", function(error, data){

		$("#search_button").click(function() {
			//console.log("hello")
			var input = $("#search").val()

			var data_length = data.length;

			for (var i = 0; i < data.length; i++) {
				if (data[i].gene == input) {
					if (data[i].file) {
						alert(data[i].file)
						break;
					}
					else {
						alert("Sorry " + input + " is in the dataset but is not super.")
						break;
					}
				}
			}

		})



	})

};

d3.selection.prototype.moveToFront = function() {
  return this.each(function(){
    this.parentNode.appendChild(this);
  });
};


var linegraph_brush, linegraph_xdomain, linegraph_ydomain;

function update_linegraph(file) {

	$(".tablesorter").trigger("update");

	//reset the clear button 
	d3.select(".clear_button_linegraph").remove();
	
	line = d3.svg.line()
		.x(function(d) {return linegraph_x(d.rank)})
		.y(function(d) {return linegraph_y(d.val)});

	//data
	d3.csv("/Documents/Bradner_work/hockey-sticks/" + file, function(error, data) {
		data.forEach(function(d) {
			d.super = +d.IS_SUPER;
			d.val = +d.SIGNAL;
			d.rank = +d.RANK;			
		});

		//console.log(data.length)

		linegraph_x.domain([0, data.length]);

		linegraph_y.domain(d3.extent(data, function(d) {
			return d.val;
		}));

		var datapoints = linegraph.selectAll(".dot")
			.data(data, function(d) {return d.rank});

		datapoints.exit().remove();


		//add rows to table
		function addRow(tableID, d) {

			var table = document.getElementById(tableID);

			//check if the dot has been added before. 
			var table_row = $("#table tbody tr");

			var rank_array = [];

			table_row.each(function(index, element) {
				var rank_row = $(this).find("td").eq(2).text();
				rank_array.push(rank_row);
			});

			//if not, add it
			if ($.inArray(d.rank, rank_array) == -1) {
				// console.log(table)

				var rowCount = table.rows.length;
				var row = table.insertRow(rowCount);

				var cell0 = row.insertCell(0);
				var element1 = document.createElement("input");
				element1.type = "checkbox";
				element1.name = "check";
				cell0.appendChild(element1);

				var celladd = row.insertCell(1)
				celladd.innerHTML = output_name[0];

				var cell1 = row.insertCell(2);
				var split = d.REGION_ID.split("_");
				if (split[3] == "lociStitched") {
					cell1.innerHTML = split[0];
				}
				else {
					cell1.innerHTML = split[3];
				}

				// console.log(d.peak)

				var cell2 = row.insertCell(3);
				cell2.innerHTML = d.RANK;

				var cell3 = row.insertCell(4);
				if (d.super == 1) {
					cell3.innerHTML = "Yes";
					cell3.style.backgroundColor = "#b30000";
				}
				else {
					cell3.innerHTML = "No";
				}

				var cell4 = row.insertCell(5);
				cell4.innerHTML = d.CHROM;

				var cell5 = row.insertCell(6);
				cell5.innerHTML = d.START;

				var cell6 = row.insertCell(7);
				cell6.innerHTML = d.STOP;

				var cell7 = row.insertCell(8);
				cell7.innerHTML = d.PROXIMAL_GENES;

				var cell8 = row.insertCell(9);
				cell8.innerHTML = d.PROXIMAL_FUNCTION;

				//highlight rows on the table on mouseover
				var rows = d3.selectAll("tbody tr")
					.on("mouseover", function() {
						d3.select(this).style("background-color", "#354299");

						var row_rank = $(this).find("td").eq(3).text();

						var here = d3.select("[rank='" + row_rank + "']")
							.moveToFront()
							.attr("r", 6.5)
							.attr("stroke", "black")
							.attr("fill", "#354299")
							.attr("stroke-width", "2px");

					})
					.on("mouseout", function() {
						d3.select(this).style("background-color", null);

						var row_rank = $(this).find("td").eq(3).text();

						var row_super = $(this).find("td").eq(4).text();

						if (row_super == "No") {
							var superness = 0;
						}
						else {
							var superness = 1;
						}

						var here = d3.select("[rank='" + row_rank + "']")
							.attr("stroke", null)
							.attr("fill", linegraph_color(superness))
							.attr("r", 4.5);
					});

				//sorting the table
				$(function() {
					$('#table').tablesorter({
						headers: {
							0: {
								sorter: false
							},
							1: {
								sorter: false
							},
							2: {
								sorter: false
							},
							5: {
								sorter: false
							},
							6: {
								sorter: false
							},
							7: {
								sorter: false
							},
							8: {
								sorter: false
							},
							9: {
								sorter: false
							}
						}
					})
				})
			};  
		}

		linegraph.select(".x.axis")
			.transition()
			.duration(1000)
			.call(linegraph_xAxis);

		linegraph.select(".y.axis")
			.transition()
			.duration(1000)
			.call(linegraph_yAxis);


		var delay = 700, clicks = 0, timer = null;

		datapoints
			.enter()
			.append("circle")
			.attr("class", "dot")
			.attr("rank", function(d) {
				return d.rank;
			})
			.attr("super", function(d) {
				return d.super;
			})
			.attr("r", 4.5)
			.attr("opacity", "0.9")
			.attr("clip-path", "url(#linegraph_clip)")
			.on("mouseover", function(d) {

				d3.select(this)
				.moveToFront()
				.attr("r", 6.5)
				.attr("stroke", "black")
				.attr("stroke-width", "2px");

				if (!d.PROXIMAL_GENES) {
					var nearby = "None";
				}
				else {
					var nearby = d.PROXIMAL_GENES;
				}
				//console.log(d)

				if (d[""]) {
					var gene = d[""];
				}
				else if(!d.top) {
					var gene = "None";
				}
				else {
					var gene = d.top;
				}

				// console.log(d[""])

				var split = d.REGION_ID.split("_");
				//console.log(d)
				if (split[3] == "lociStitched") {
					var peak = split[0];
				}
				else {
					var peak = split[3];
				}

				graph_tip.html("<strong>Peak Number: </strong>" + peak + "<br><strong>Top Gene: </strong>" + gene + "<br><strong>Chromosome: </strong>" + d.CHROM + "<br><strong>Nearby genes: </strong>" + nearby);
				graph_tip.show(d);
			})
			.on("mouseout", function(d) {
				d3.select(this)
				.attr("r", 4.5)
				.attr("stroke", null)

				graph_tip.hide(d);
			})
			.on("click", function(d) {
				clicks++;

				//console.log(d)

				if (clicks === 1 & d.super == 1) {

					//console.log(file)
//output_name[0] + "plots/" +
					var pdf_name = "/Documents/Bradner_work/hockey-sticks/" + output_name[0] + "_plots/" + "SE_plots_" + output_name[0] + "_" + d.REGION_ID + ".pdf";

					//console.log(output_name)

					//console.log(pdf_name)
					//U87_H3K27AC_plots/SE_plots_U87_H3K27AC_1_MACS_peak_6695_lociStitched.pdf
					//U87_H3K27AC_plots/SE_plots_U87_H3K27AC6_MACS_peak_4435_lociStitched.pdf

					timer = setTimeout(function() {
						clicks = 0;
						//console.log("happen")
						screenshotPreview(pdf_name);

						// d3.selectAll(".pdf_image", function(d) {
						// 	console.log("here")
						// 	this.attr("transform", "translate(" + d.rank + "0)");
						// })
						
					}, delay);
				}

				else if (clicks === 1) {
					timer = setTimeout(function() {
						clicks = 0;
					}, delay);
				}

				else {

					clearTimeout(timer);
					clicks = 0
					addRow("tbody", d);
				}
			});

			d3.select("#linegraph")
				.on("click", function() {
					d3.select(".pdf_image").remove();
				})

		// console.log(data)

		datapoints
			.attr("cx", linegraph_x(1))
			.attr("cy", linegraph_y(1))
			.transition()
			.duration(1000)
			.attr("cx", function(d) {
				return linegraph_x(d.rank);
			})
			.attr("cy", function(d) {
				return linegraph_y(d.val);
			})
			.attr("fill", function(d) {
				return linegraph_color(d.super);
			});


		//set up the color legend
		var legend = linegraph.append("g")
			.attr("class", "legend")
			.attr("width", 200)
			.attr("height", 200)
			.selectAll("g")
			.data(linegraph_color.domain().slice().reverse())
			.enter()
			.append("g")
			.attr("transform", function(d,i) {
				return "translate(-10," + -1.5*(i*20) + ")";
			});

		legend.append("rect")
			.attr("width", 18)
			.attr("height", 18)
			.style("fill", linegraph_color);

		legend.append("text")
			.attr("x", 25)
			.attr("y", 0)
			.attr("dy", "1.15em")
			.text(function(d) {
				if (d == 1) {
					return "Super-enhancer"
				}
				else {
					return "Typical Enhancer"
				}
			});

		linegraph.append("text")
			.attr("x", 25)
			.attr("y", -40)
			.text("Legend")
			.style("font-weight", "bold")
			.style("font-size", "14px");

		//set up the brush
		linegraph_brush = d3.svg.brush()
			.x(linegraph_x)
			.y(linegraph_y)
			.on("brushend", linegraph_brushend);

		var max_x = d3.max(data, function(d) {
			return d.rank;
		});
		var max_y = d3.max(data, function(d) {
			return d.val;
		});

		//console.log(max_y)

		linegraph_xdomain = [0, max_x];
		linegraph_ydomain = [0, max_y];

		linegraph_x.domain(linegraph_xdomain);
		linegraph_y.domain(linegraph_ydomain);

		d3.select(".linegraph_brush").call(linegraph_brush);

		//updating the brush based on selected data

		var linegraph_clear_button;

		function linegraph_brushend() {

			var linegraph_x_domain = [linegraph_brush.extent()[0][0], linegraph_brush.extent()[1][0]]
			var linegraph_y_domain = [linegraph_brush.extent()[0][1], linegraph_brush.extent()[1][1]]

			// console.log(linegraph_x_domain)
			// console.log(linegraph_y_domain)

			// equal domain ends means click on graph
			// coerce dates to numbers to check equality
			if (+linegraph_x_domain[0] == +linegraph_x_domain[1] && +linegraph_y_domain[0] == +linegraph_y_domain[1]) {
				return;
			}

			get_button = d3.select(".linegraph_clear_button");

			if (get_button.empty() === true)
			{
				linegraph_clear_button = linegraph.append("g")
					.attr("transform", "translate(" + (bb_linegraph.w - 100) + "," + (bb_linegraph.h - 320) + ")")
					.attr("class", "linegraph_clear_button");

				linegraph_clear_button.append("rect")
					.attr("width", 102)
					.attr("height", 20)
					.attr("y", -17)
					.attr("x", -4)
					.attr("rx", "10px")
					.attr("ry", "10px")
					.style("fill", "#9f9f9f");

				linegraph_clear_button
					.append('text')
					.attr("y", -1)
					.attr("x", 10)
					.text("Clear Zoom")
					.style("fill", "black");
			}

			linegraph_x.domain(linegraph_x_domain);
			linegraph_y.domain(linegraph_y_domain);

			linegraph_transition();

			d3.select(".linegraph_brush").call(linegraph_brush.clear());

			// add the on click events for the button
			linegraph_clear_button.on('click', function ()
			{
			    // reset everything
				linegraph_x.domain(linegraph_xdomain);
				linegraph_y.domain(linegraph_ydomain);

				linegraph_transition();

				linegraph_clear_button.remove();
			});

			function linegraph_transition() {
				
				linegraph.select(".x.axis")
					.transition()
					.duration(1000)
					.call(linegraph_xAxis);

				linegraph.select(".y.axis")
					.transition()
					.duration(1000)
					.call(linegraph_yAxis);

				linegraph.selectAll("circle")
					.transition()
					.duration(1000)
					.attr("cx", function(d) {
						return linegraph_x(d.rank);
					})
					.attr("cy", function(d) {
						return linegraph_y(d.val);
					});
			}
		}
	});
}

function deleteRow(tableID) {
	var table = document.getElementById(tableID);
	var rowCount = table.rows.length;

	for (var i = 0; i < rowCount; i++) {
		var row = table.rows[i];
		var checkbox = row.cells[0].childNodes[0];

		if (null!= checkbox && checkbox.checked == true) {
			table.deleteRow(i);
			rowCount--;
			i--;
		}
	}

	$(".tablesorter").trigger("update");
}	

function clearTable(tableID) {
	var table = document.getElementById(tableID);
	var rowCount = table.rows.length;

	for (var i = 0; i < rowCount; i++) {
		table.deleteRow(i);
		rowCount--;
		i--;
	}
}

//export to excel file on click of export button
$(document).ready(function () {

    function exportTableToCSV($table, filename) {

    	// console.log($table)

        var $rows = $table.find('tr:has(td)');
        // console.log($rows)

            // Temporary delimiter characters unlikely to be typed by keyboard
            // This is to avoid accidentally splitting the actual contents
            tmpColDelim = String.fromCharCode(11); // vertical tab character
            tmpRowDelim = String.fromCharCode(0); // null character

            // actual delimiter characters for CSV format
            colDelim = ',';
            rowDelim = '"\r\n"';

            // Grab text from table into CSV formatted string
            csv = '"' + $rows.map(function (i, row) {
                var $row = $(row),
                    $cols = $row.find('td');

                return $cols.map(function (j, col) {
                    var $col = $(col),
                        text = $col.text();

                    return text.replace('"', '""'); // escape double quotes

                }).get().join(tmpColDelim);

            }).get().join(tmpRowDelim)
                .split(tmpRowDelim).join(rowDelim)
                .split(tmpColDelim).join(colDelim) + '"';

            // Data URI
            csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

        $(this)
            .attr({
            'download': filename,
                'href': csvData,
                'target': '_blank'
        });
    }

    // This must be a hyperlink
    $("#export_button").on('click', function (event) {
        // CSV
        exportTableToCSV.apply(this, [$('#table_div>table'), 'export.csv']);

    });	
});


