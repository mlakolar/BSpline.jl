module AquaTest

using Test
using BSpline
using Aqua

@testset "Aqua" begin
  Aqua.test_all(BSpline)
end

end
