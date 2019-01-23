########################################################################################
#
# erode.jl
#
########################################################################################


using DelimitedFiles
using SparseArrays


### data types


mutable struct Cell                     # a location on surface
    x::Float64                          # center coordinate
    y::Float64
    z::Float64                          # elevation (subject to change)
	A::Float64 							# drainage area
	grad::Float64 						# mean gradient magnitude
	lake::Bool 							# is part of a lake
    neighbors::Array{Int64, 1}          # list of connecting cells
end


mutable struct Connection               # cell connections
    cell1::Int64                        # connecting cell index numbers
    cell2::Int64
    deltaZ::Float64                     # cell_1 elevation - cell_2 elevation
    w::Float64                          # connection weight
end


mutable struct Params
    d::Float64                          # grid spacing (uniform in x & y)
	zMin::Float64 						# minimum elevation; extracted from edges
    scaleMin::Float64                   # scaling factors used to weight impact of drainage area on erosion rate
	scaleMax::Float64
	scaleF::Float64
	dzMax::Float64 						# maximum elevation changer per iteration
	iterMax::Int64 						# number of iterations; each iteration = step that removes dzMax
    maxSinks::Int64                     # maximum number of internal drainage cells allowed in domain
    iterMaxFill::Int64                  # maximum number of internal drainage fill attempts, per outer iteration
	gradLake::Float64 					# maximum gradient magnitude threshold to include cell within a lake feature
    elevLake::Float64                   # lake elevation elevation threshold (alternative means for defining lake)
	lakeCutoff::Float64 				# threshold distance for including/excluding lake points in lake group
end


### input and output functions


function ReadSurface()
    # read surface properties and populate array of cell structs
    cell = Cell[]
    data = readdlm("surface_input.csv", ',', header=true)
    xArray = data[1][:, 1]                          # determine extent of x-y grid
    yArray = data[1][:, 2]
    nx = length(collect(Set(xArray)))
    ny = length(collect(Set(xArray)))    
    for i = 1:size(data[1], 1)
        x = Float64(data[1][i, 1])                  # must be sorted by x (inner) and then y (outer)
        y = Float64(data[1][i, 2])
        z = Float64(data[1][i, 3])
		A = 0.0 									# drainage area placeholder
 		grad = 0.0 									# gradient magnitude placeholder
		lake = false 									# default assumption
        push!(cell, Cell(x, y, z, A, grad, lake, Int64[]))      # connecting cells list to be populated elsewhere
    end
    println("Read surface.")
    return cell, nx, ny
end


function ReadParams()::Params
    # read numerical model parameters
    data = readdlm("params.txt", '\t', header=false)
    d = Float64(data[1, 2])
    scaleMin = Float64(data[2, 2])
    scaleMax = Float64(data[3, 2])	
    dzMax = Float64(data[4, 2])	
	iterMax = Int64(data[5, 2])
	maxSinks = Int64(data[6, 2])
	iterMaxFill = Int64(data[7, 2])
	gradLake = Float64(data[8, 2])
	elevLake = Float64(data[9, 2])    
	lakeCutoff = Float64(data[10, 2])
	scaleF = 0. 					# placeholders; defined dynamically within code
	zMin = 0.
    params = Params(d, zMin, scaleMin, scaleMax, scaleF, dzMax, iterMax, maxSinks, iterMaxFill, gradLake, elevLake, lakeCutoff)
    println("Read model parameters.")
    return params
end


function WriteCells(cell::Array{Cell, 1}, fileName::String)
    # summarize cell properties to file
    csvfile = open(fileName,"w")
    line_out = "x" * "," * "y" * "," * "z" * "," * "drainage" * "," * "gradient"
    println(csvfile, line_out)
    for (i, ce) in enumerate(cell)
		line_out = string(ce.x) * "," * string(ce.y) * "," * string(ce.z) * "," * string(ce.A) * "," * string(ce.grad)
		println(csvfile, line_out)
    end
    close(csvfile)
    println("Wrote " * fileName * ".")
end


function WriteLakePerim(lakeContainer)
	# write lake perimeter points to file (for use in flow & transport model input, etc.)
    csvfile = open("LakePerimeters.csv","w")
    line_out = "x" * "," * "y" * "," * "lake"
    println(csvfile, line_out)
    for (i, lakePerim) in enumerate(lakeContainer)
		for j = 1:length(lakePerim[1])
			line_out = string(lakePerim[1][j]) * "," * string(lakePerim[2][j]) * "," * string(i)
			println(csvfile, line_out)
		end
    end
    close(csvfile)
    println("Wrote lake perimeters file.")	
end


### cell connectivity and area + erosion modeling functions


function ElevMin(cell::Array{Cell, 1})::Float64
    # find minimum starting elevation; this will fix the lowest allowed elevation anywhere
    zMin = 1e+10
    for ce in cell
        zMin = minimum([zMin, ce.z])
    end
    return zMin
end


function Connection(cell::Array{Cell, 1}, nx::Int64, ny::Int64)
    # populate cell connection arrays
    connect = Connection[]
	deltaZ = 0.0 					# placeholder values
	w = 0.0
    for j = 1:ny, i = 1:nx-1                # x-direction connections
        cell1 = i + (j-1)*nx
        cell2 = cell1 + 1
        push!(cell[cell1].neighbors, cell2)
        push!(cell[cell2].neighbors, cell1)
        push!(connect, Connection(cell1, cell2, deltaZ, w))
    end
    for i = 1:nx, j = 1:ny-1                # y-direction connections
        cell1 = i + (j-1)*nx
        cell2 = i + j*nx
        push!(cell[cell1].neighbors, cell2)
        push!(cell[cell2].neighbors, cell1)
        push!(connect, Connection(cell1, cell2, deltaZ, w))        
    end        
    return cell, connect
end


function Outflow(cell::Array{Cell, 1})
    # calculate net outflows from cells
    cellOut = Float64[]
    sinks = Int64[]
	for (i, ce) in enumerate(cell)
		outflowTot = 0.0
		for nghbr in ce.neighbors
			outflowTot += (ce.z+0.0-cell[nghbr].z) * ((ce.z+0.0-cell[nghbr].z)>0.0)
		end
        push!(cellOut, outflowTot)
        if outflowTot==0.0      # ce currently is an (unallowed) internal sink
            push!(sinks, i)
        end
	end
    return cellOut, sinks
end


function Fill(cell::Array{Cell, 1}, params::Params)
    # fill in all internal sinks
    cellOut, sinks = Outflow(cell)
    iCounter = 0
    while (length(sinks) > params.maxSinks)  && (iCounter < params.iterMaxFill)
		iCounter += 1
        for s in sinks
            zTot = 0.0
            for nghbr in cell[s].neighbors
                zTot += cell[nghbr].z
            end
			cell[s].z = zTot/length(cell[s].neighbors)
        end
        cellOut, sinks = Outflow(cell)
    end
    return cell, cellOut, length(sinks)
end


function Weights(cell::Array{Cell, 1}, cellOut::Array{Float64, 1}, connect::Array{Connection, 1})::Array{Connection, 1}
	# allocate outflow drainage area weights across connections
    for cn in connect
        cn.deltaZ = cell[cn.cell1].z - cell[cn.cell2].z
        if cn.deltaZ > 0.0
            cn.w = cn.deltaZ/cellOut[cn.cell1]
        elseif cn.deltaZ < 0.0
            cn.w = abs(cn.deltaZ)/cellOut[cn.cell2]        
        else
            cn.w = 0.0
        end
    end
	return connect
end


function Gradient(cell::Array{Cell, 1}, params::Params)::Array{Cell, 1}
	# magnitude of the land surface gradient, per cell
	for ce in cell
		xSum = 0.0
		ySum = 0.0
		numX = 0
		numY = 0
		for nghbr in ce.neighbors
			xSum += (ce.y == cell[nghbr].y) * (ce.z - cell[nghbr].z) * sign(ce.x - cell[nghbr].x)
			numX += (ce.y == cell[nghbr].y) * 1
			ySum += (ce.x == cell[nghbr].x) * (ce.z - cell[nghbr].z) * sign(ce.y - cell[nghbr].y)
			numY += (ce.x == cell[nghbr].x) * 1			
		end
		xGrad = xSum/(numX * params.d)
		yGrad = ySum/(numY * params.d)
		ce.grad	= sqrt(xGrad^2 + yGrad^2)	
	end
	return cell
end


function AssembleMatrix(cell::Array{Cell, 1}, connect::Array{Connection, 1})
    # fill out the LHS of the equation matrix
    row_index = Int64[]                     # indexing system for sparse matrix
    col_index = Int64[]
    data = Float64[]
    for i = 1:length(cell)              # set diagonal term
        push!(row_index, i)         
        push!(col_index, i)
        push!(data, 1.0) 
    end
    # off-diagonal terms from across each connection; matrix will be asymmetric
    for cn in connect
        if cn.deltaZ > 0.0      # flow into cell2
            push!(row_index, cn.cell2)
            push!(col_index, cn.cell1)
            push!(data, -cn.w)
        elseif cn.deltaZ < 0.0      # flow into cell1
            push!(row_index, cn.cell1)
            push!(col_index, cn.cell2)
            push!(data, -cn.w)
		end
    end
    b = zeros(length(cell)).+ 1.0
    return data, row_index, col_index, b
end


function Chisel(cell::Array{Cell, 1}, params::Params)::Array{Cell, 1}
    # calculate elevation changes at each cell
    dz0 = Float64[]
    for ce in cell
        push!(dz0, -ce.grad * (ce.A^params.scaleF))
    end
    dzn = maximum(abs.(dz0))
    dz = params.dzMax * dz0/dzn
    for i = 1:length(cell)
        cell[i].z += dz[i] * ((cell[i].z+dz[i])>params.zMin)
    end
    return cell
end
    

function LakeCells(cell::Array{Cell, 1}, params::Params)
	# delineate subset of cells within lake(s)
	lakeCell = Cell[]
	for ce in cell
        if (ce.grad <= params.gradLake) || (ce.z <= params.elevLake)
			push!(lakeCell, ce)
			ce.lake = true
		end
	end
	return lakeCell, cell
end


function GroupLakes(lakeCell::Array{Cell, 1}, params::Params)
	# group lake cells into clusters; return arrays of points
	lakeGroups = []
	push!(lakeGroups, Int64[])
	push!(lakeGroups[end], 1)
	for i = 2:length(lakeCell) 								# for each point
		match = false
		j = 1
		while (j<=length(lakeGroups)) && (match==false) 					# for each lake group
			k = 1
			while (k<=length(lakeGroups[j])) && (match==false) 			# for each lake group member
				lakePt = lakeGroups[j][k] 						# lake cell index number corresponding to j & k
				dist = sqrt((lakeCell[i].x-lakeCell[lakePt].x)^2 + (lakeCell[i].y-lakeCell[lakePt].y)^2)
				if dist <= params.lakeCutoff
					push!(lakeGroups[j], i) 		# add lake index to current group
					match = true
					break
				end
				k += 1
			end
			j += 1
		end
		# fallen thru filter - point belongs to no lake; create new one
		if match==false
			push!(lakeGroups, Int64[])
			push!(lakeGroups[end], i)
		end		
	end
	return lakeGroups
end
	

function Perimeter(lakeCell::Array{Cell, 1}, lakeGroup::Array{Int64, 1})

	# perform simple bubble sort of lakeGroup list, by x-coordinate of each lake cell
	if length(lakeGroup)>1
		ordered=false 			# default assumption
		while ordered == false
			ordered = true
			for i = 2:length(lakeGroup)
				if lakeCell[lakeGroup[i]].x < lakeCell[lakeGroup[i-1]].x
					# these two locations are listed in decreasing order
					ordered = false
					a = lakeGroup[i]
					lakeGroup[i] = lakeGroup[i-1]
					lakeGroup[i-1] = a
				end	
			end
		end
	end
	
	# determine perimeter
	xp = Float64[]
	yp = Float64[]
	x = lakeCell[lakeGroup[1]].x
	yMin = lakeCell[lakeGroup[1]].y
	yMax = yMin
	numCols = 1
	for i = 2:length(lakeGroup)
		if lakeCell[lakeGroup[i]].x != lakeCell[lakeGroup[i-1]].x
			# starting a new column; add old extrema to xp and yp arrays
			insert!(xp, numCols, lakeCell[lakeGroup[i-1]].x) 	# bottom location
			insert!(yp, numCols, yMin)
			insert!(xp, numCols+1, lakeCell[lakeGroup[i-1]].x) 	# top location
			insert!(yp, numCols+1, yMax)			
			numCols += 1
			yMin = 1e+10
			yMax = -1e+10
		end
		yMin = minimum([yMin, lakeCell[lakeGroup[i]].y])
		yMax = maximum([yMax, lakeCell[lakeGroup[i]].y])		
	end
	# add extrema from last column to arrays
	insert!(xp, numCols, lakeCell[lakeGroup[end]].x) 	# bottom location
	insert!(yp, numCols, yMin)
	insert!(xp, numCols+1, lakeCell[lakeGroup[end]].x) 	# top location
	insert!(yp, numCols+1, yMax)
	
	return xp, yp
end

	
### main script


function Erode()

    # read and process input files
    cell, nx, ny = ReadSurface()                            # read surface properties
    params = ReadParams()                                   # read model parameters
    
    # process initial surface definition
    params.zMin = ElevMin(cell)
    cell, connect = Connection(cell, nx, ny)                # tabulate cell connections
	println("Processed cell network.")

    # iteratively fill in internal drainages       
    cell, cellOut, numSinks = Fill(cell, params)
    
	sinks = []
    for i = 1:params.iterMax
    
        println("\titeration " * string(i) * "; no. of sinks = " * string(numSinks))    
        params.scaleF = params.scaleMin + (i/params.iterMax) * (params.scaleMax - params.scaleMin)
        
        # assign area distribution weights, by connection
        connect = Weights(cell, cellOut, connect)        
        
        # assemble matrix of flow balance equations
        data, row_index, col_index, b = AssembleMatrix(cell, connect)
        A = sparse(row_index, col_index, data, length(cell), length(cell))
        
        # solve system of equations
        area = \(A, b)
            
        # update surface parameters: drainage area and local gradient magnitude
        for i = 1:length(cell)
            cell[i].A = area[i] * params.d^2
        end
        cell = Gradient(cell, params)  

        # simulate localized preferential erosion
        cell = Chisel(cell, params)
        
        # iteratively fill in internal drainages
        cell, cellOut, numSinks = Fill(cell, params)
    
    end    

	# process (implied) lakes
	lakeCell, cell = LakeCells(cell, params) 		# note locations of internal lake cells (by low gradient)
	lakeGroups = GroupLakes(lakeCell, params) 		# array of lake index groups, corresponding to lakes
	lakeContainer = [] 								# determine lake outline perimeters
	for (i, lakes) in enumerate(lakeGroups)
		println("Processing lake no. " * string(i) * "; cells = " * string(length(lakes)))
		xp, yp = Perimeter(lakeCell, lakes)
		push!(lakeContainer, [xp, yp])
	end
	WriteLakePerim(lakeContainer)
	WriteCells(lakeCell, "Lakes.csv")	# write lake cell sets to separate file (for plotting, etc.)



    # delineate streams; need to figure out if cell that is being crossed has its lake flag checked
    
    # write elevation and drainage area results
    WriteCells(cell, "FinalState.csv")           # write cell output file
    
    println("Finished.")
    
end

### run script
Erode()