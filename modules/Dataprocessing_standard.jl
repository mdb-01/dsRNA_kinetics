""" 
    slope_estimation_standard(data::DataFrame)

    Calculates the slope of converting intensity from dotblot measurement to concentration of dsRNA (ng/mL)
    from standard signal
"""

function slope_estimation_standard(data::DataFrame,  sigma::Float64, lambda::Float64=1e-5)
    # data includes a column with the known concentrations and a column with the respective volume 

    #initialise storage df with slope and uncertainty
    slope_error_intersect = zeros(1,3) #column 1 = slope, column 2 = uncertainty, column 3 = intersection

    n_points = size(data, 1) #number of datapoints
    
    X = [ones(n_points) data.concentration]
    phi = sigma.^2*I(n_points)
    Y = data.volume

    #regularisation to improve stability
    new_matrix = transpose(X) * inv(phi) * X #+ lambda * I(size(X, 2))

    # Check condition number to ensure stability
    if cond(new_matrix) > 1e10  # If condition number is large, matrix is ill-conditioned
        println("Matrix is nearly singular, using pseudo-inverse")
        B = pinv(new_matrix) * transpose(X) * inv(phi) * Y  # Use pseudo-inverse instead of inverse
    else
        B = inv(new_matrix) * transpose(X) * inv(phi) * Y   # Normal inversion
    end

       
    #B = inv(transpose(X)*inv(phi)*X)*transpose(X)*inv(phi)*Y        # find B0 and B1 -> Y = B0 + B1*X
                                                                  
    slope_error_intersect[1,1] = B[2]                                # store B1 (eg slope)
    slope_error_intersect[1,3] = B[1]                                # store intersection point
        
    cov = inv(transpose(X)*inv(phi)*X)                              # calculate covariance of b #[cov]=
    concentrations = data.concentration
    variance = [[x,1]'*cov*[x,1] for x in concentrations]           # calculate the variance 
    variance_average = sum(variance)                                # average across time points
    slope_error_intersect[1,2] = sqrt(variance_average)             # store the st dev of the replicate



    return slope_error_intersect, variance, cov
end     


function initial_fit(data::DataFrame)
    n_points = size(data, 1)
    X = [ones(n_points) data.concentration]
    Y = data.volume
    
    B = inv(transpose(X) * X) * transpose(X) * Y  # Least squares estimation
    fitted_values = X * B
    residuals = Y - fitted_values
    return residuals
end

# Then, use the residuals to estimate sigma
function estimate_sigma(residuals::Vector)
    return std(residuals)  # Standard deviation of residuals
end

function initial_fit_relative(data::DataFrame)
    n_points = size(data, 1)
    X = [ones(n_points) data.concentration]
    Y = data.volume
    
    # Least squares estimation
    B = inv(transpose(X) * X) * transpose(X) * Y  
    fitted_values = X * B

    # Calculate relative residuals
    residuals_relative = (Y - fitted_values) ./ fitted_values  # Relative residuals as a percentage
    return residuals_relative
end
