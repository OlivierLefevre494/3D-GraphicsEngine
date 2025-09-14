mainframe = document.getElementById('mainframe')
box = document.getElementById('box')

class Vector {
    constructor(...components) {
        this.components = components
    }

    add({components}){ // works
        (components.map((component, index) => this.components[index] += component))
        return this
    }

    static sum(vector1, vector2) {
        var tmp = new Vector()
        vector1.components.map((component, index) => {tmp.components.push(vector2.components[index] + component)})
        return tmp
    }

    iszero(){ // works even with float 0.0
        var iszero = true
        this.components.map((component, index) => {if (component!=0) {iszero=false}})
        return iszero
    }

    dot(vector) { // works
        var sum = 0
        this.components.map((component, index) => {sum += ((component * vector.components[index]))})
        return sum
    }

    subtract(components) { // works
        components.components.map((component, index) => {this.components[index] -= component})
        return this
    }

    static difference(vector1, vector2) { // from vector 1 to vector 2 A->B works
        var tmp = new Vector()
        tmp.components = []
        vector1.components.map((component, index) => {tmp.components.push(vector2.components[index] - component)})
        return tmp
    }

    static islinear(vector1, vector2) { // works
        // this is a static function to check if two vector's are linear combinations of each other, same up to scaling
        if ((vector1.components.length != vector2.components.length) | ((vector1.components.length * vector2.components.length) == 0)) {
            console.error('Either your vectors have not matching number of components or one of them is 0.')
        } else {
            var coef = 0
            vector1.components.map((component, index ) => { if ((component!=0) & (vector2.components[index]!=0)) {coef = component/vector2.components[index]}})
            var linear = true
            vector1.components.map((component, index) => { if ((vector2.components[index] * coef)!=component) {linear = false} })
            return linear
        }
    }

    midpoint() {
        return null // add function to get midpoint this will probably be useful...
    }

    multiply(scalar) {
        var tmp = new Vector()
        tmp.components = []
        this.components.map((component, index) => {tmp.components.push(component * scalar)})
        return tmp
    }

    static isperpendicular(vector1, vector2) { // works
        return (vector1.dot(vector2) == 0)
    }

    magnitude() { // works
        var sum = 0
        this.components.map((component, index) => {sum += Math.pow(component, 2)})
        return Math.pow(sum, 1/2)
    }

    proj(vector) { // given an instance, and another vector, project this vector onto that one
        var tmp = new Vector()
        tmp.components = JSON.parse(JSON.stringify(vector.components))
        return [(tmp.multiply(tmp.dot(this) / (Math.pow(vector.magnitude(),2)))),(tmp.dot(this) / ((Math.pow(vector.magnitude(),2))))]  // a . b = cos(ß) * ||a|| * ||b|| = ? * ||b|| / ||b||^2 = ?/||b|| * b = ?
    }

    static getperp(vector) {
        // if (orientation == 'r') {
        var tmp = new Vector(-vector.components[1], vector.components[0])
        // } else if (orientation == 'l') {
        //     var tmp = new Vector(vector.components[1], -vector.components[0])
        // }
        return tmp
    }

    static getangle(vec1, vec2, orientation) {
        var perp = Vector.getperp(vec1)
        var xfac
        var x
        var tmp = vec2.proj(perp)
        x = tmp[0]
        xfac = tmp[1]
        var yfac
        var y
        var tmp = vec2.proj(vec1)
        y = tmp[0]
        yfac = tmp[1]
        var cosangle = Math.acos((Math.sign(xfac)*(x.magnitude()))/vec2.magnitude())/Math.PI
        var sinangle = 0.5-Math.asin((Math.sign(yfac)*(y.magnitude()))/vec2.magnitude())/Math.PI
        var totalangle
        console.log(sinangle)
        if ((sinangle<=1/2)) {
                totalangle = 0.5+cosangle
            } else {
                totalangle = 1/2-cosangle
                if (totalangle<0) {
                    totalangle = 2+totalangle
                }
            }
        if (orientation=='r') {
            totalangle = 2-totalangle
        }
        return totalangle
    }

}

class Vertex {

    constructor(...components) {
        this.components = components
        this.screencomponents = []
        this.neighbours = []
    }

    addneighbour(...neighbours) {
        neighbours.map((neighbour, index) => {this.neighbours.push(neighbour)})
    }
    
    
}

class Box {
    constructor (...vertices) {
        this.vertices = vertices
    }
}

class Square {

    constructor (...vertices) { // works

        this.vertices = vertices
        if (this.vertices.length != 4) {
            throw new Error('Not enough vertices to make a square!')
        }
        // Make sure opposite vectors are linear
        this.segments = []
        this.vertices.map((vertex, index) => {this.segments.push(Vector.difference(vertex, this.vertices[((index + 1)%4)]))})
        // I think I only need to check if they are perpendicular but this code is nice
        this.segments.map((segment, index) => {if (!Vector.islinear(segment, this.segments[((index + 2)%4)])) {throw new Error('Not squarey')}})
        // checking perpendicularity
        this.segments.map((segment, index) => {if (!Vector.isperpendicular(segment, this.segments[((index + 1) % 4)])) {throw new Error('Not perpendicalified rip bozo')}})
    }

}

class World {
    constructor(dim=3) {
        this.dimensions = dim
        this.objects = []
    }

    addobject(object) {
        this.objects.push(object)
    }

    project(observer){
        var previousprojections = document.getElementsByClassName('projection')
        for (let i=0; i<previousprojections.length; i++) {
            previousprojections[0].remove()
        }
        this.objects.map((object, index) => {object.display(observer)})
    }


}


/* Creating a projection should work, possible problems : sizing. Should work apart from that. Should check if same point or if diff < something*/


function createProjection(point1, point2, pixelsx=innerWidth, pixelsy=innerWidth) {
    console.log("PROJECTIONS")
    console.log("POINTS 1 AND TWO")
    console.log(point1)
    console.log(point2)
    new_disp = document.createElement('div')
    new_disp.classList = 'projection'
    new_disp.style.height = '1px'
    new_disp.style['background'] = 'black'
    new_disp.style['position'] = 'absolute'

    new_disp.style['z-index'] = 1
    var diff = Vector.difference( point1,point2) // the first component should be the x coords
    var angle = (Math.atan(diff.components[1]/diff.components[0]))
    new_disp.style.transformOrigin = 'bottom left'
    new_disp.style.rotate = String(-angle % 2) + 'rad' // create a point with an angle!
    diff.components[0] = diff.components[0] * pixelsx // this will get pixels in x dir for mag
    diff.components[1] = diff.components[1] * pixelsy // this will get the pixels in y dir 
    var mag = diff.magnitude() // this will be in pixels now, cool! // - Math.abs((mag/2)*(1-Math.cos(angle))) +  + (Math.sin(angle) * (mag/2))
    new_disp.style['left'] = String((point1.components[0] * pixelsx)) + 'px'
    new_disp.style['bottom'] = String((point1.components[1] * pixelsy)) + 'px'
    new_disp.style['width'] = String(mag) + 'px' // get magnitude of delta vec
    new_disp.class = 'projection'
    document.body.appendChild(new_disp)
    // get two points, then get angle then get length and then make line easy as pie
}

function arraysEqual(arr1, arr2) {
    if (arr1.length != arr2.length) {return false}
    for (let i=0; i<arr1.length; i++) {
        if (arr1[i]!=arr2[i]) {
            return false
        }
    }
    return true
}

class object{
    constructor(...vertices) {
        this.vertices = vertices
        this.faces = []
    }

    display(observer) {
        // console.log("DISPLAYING SHIT")
        // console.log(this.vertices)
        var vectors = this.vertices.map((vertex, index) => {return (Vector.difference(vertex,observer.position))}) // get all vectors
        // console.log("VECTORS RELATIVE TO OBSERVER POSITION")
        // console.log(observer.position)
        // console.log(vectors)
        // console.log("Vectors dot / screen distance")
        // console.log(vectors.map((vector, index) => {return vector.dot(observer.direction)}))
        // console.log(vectors.map((vector, index) => {return vector.dot(observer.direction)/Math.pow(observer.screendistance,1)}))
        // console.log("DIRECTION")
        // console.log(observer.direction)
        // console.log("Vectors multiplied")
        // console.log(vectors.map((vector, index) => {return vector.multiply(Math.pow(observer.screendistance,1)/vector.dot(observer.direction))}))
        // console.log("CORNERS")
        // console.log(observer.screen.corners)
        var points_on_screen = vectors.map((vector,index) => {return (vector.multiply(Math.pow(observer.screendistance,1)/vector.dot(observer.direction)).subtract(observer.screen.corners[2])) /* This will get points on screen, sub bottom left corner*/})
        // console.log("POINTS ON SCREEN")
        // console.log(points_on_screen)
        var screen_basis_of_points = points_on_screen.map((point, index) => {return observer.screen.basischange(point)}) // implement method in screen to get new coordinates under the basis
        this.vertices.map((vertex, index) => {this.vertices[index].screencomponents[0] = screen_basis_of_points[index].components[0] * innerWidth; this.vertices[index].screencomponents[1]=screen_basis_of_points[index].components[1]*innerWidth;})
        // console.log("SCREEN BASIS OF POINTS")
        // console.log(screen_basis_of_points)
        // now all these points have neighbours
        var created_edges  = []
        for (let i=0; i<screen_basis_of_points.length; i++) {
            var index = i
            var point = screen_basis_of_points[i]
            for (let j=0; j<this.vertices[index].neighbours.length; j++) {
                var neighbour = this.vertices[index].neighbours[j]
                var idx = 0
                this.vertices.map((value, index) => {if (value==neighbour) {idx=index}})
                var neighbour_point = screen_basis_of_points[idx]
                var alreadybuilt = false
                for (let k=0; k<created_edges.length; k++) {
                    if (arraysEqual(created_edges[k], [neighbour_point, point])) {
                        alreadybuilt = true
                    }
                }
                if (alreadybuilt) { // if already in list idk contains or whatever
                } else {
                    created_edges.push([point, neighbour_point])
                    created_edges.push([neighbour_point, point])
                    if (point.components[0]<=neighbour_point.components[0]) {
                        createProjection(point, neighbour_point)
                    } else {
                        createProjection(neighbour_point, point)
                    }
                
                }
            }

        } // so now I've projected onto the screen and I am done!
        /* Essentially here I am getting all of the coordinates of the vertices relative to the observer
        Then I get the points on screen by dividing it by the scale (trigonometry) and subtracting the bottom left corner, then with the dot product for each axis
        Then I store those in the object under like coordinates, then I can project If I want oh wait this is the project function lol. 
        */
        // I have all the points relative to screen, but now I need to get all the vertices through some sort of graph traversal in a way that makes sense
        // I.E you have points ...... need to pair them to get vectors ... which you can project on the screen. I can display on screen, I have the coordinates i 
        this.faces.map((face, index)=>{face.display()})
    }

}

class face {
    constructor(color, ...vertices) { // this should have vertices s.t a -> b-> c-> d -> a.... In an order
        this.color = color
        this.vertices = vertices
        this.screen_vertices = null
        this.orientation = ''
    }

    checkorientation() { // check what is inside the shape
        var left = 0
        var right = 0
        for (let i=0; i<this.screen_vertices.length;i++) {
            var vec1 = Vector.difference(this.screen_vertices[i], this.screen_vertices[(i+1)%this.screen_vertices.length])
            var vec2 = Vector.difference(this.screen_vertices[(i+1) % this.screen_vertices.length], this.screen_vertices[(i+2) % this.screen_vertices.length])
            var angle = Vector.getangle(vec1, vec2, 'r')
            right += angle
            left += 2-angle
        }
        if (left>right) { // right orientation is 'inside' the shape
            return 'r'
        } else if (left<right) {
            return 'l'
        } else {
            return Error("Left and right orientation equal, shape musn't close")
        }

    }

    display() {
        
        var copiedvertices = this.copy()
        this.screen_vertices = copiedvertices
        this.orientation = this.checkorientation()
        var tmpunique = [...Array(this.vertices.length).keys()]
        var i = 0
        while ((3<=tmpunique.length) & (i<10)) {
            var firstvec
            var secondvec
            var thirdvec
            var midtocenter
            var center
            console.log(copiedvertices[tmpunique[i%tmpunique.length]])
            console.log("BEGINNING PROCESS ------------------------------------- CHECKING THE TRIANGLE WITH VERTICES ")
            console.log(copiedvertices[tmpunique[i%tmpunique.length]])
            console.log(copiedvertices[tmpunique[(i + 1)%tmpunique.length]])
            console.log(copiedvertices[tmpunique[(i + 2)%tmpunique.length]])
            console.log(copiedvertices)
            var a = copiedvertices[tmpunique[i%tmpunique.length]].components
            var b = copiedvertices[tmpunique[(i + 1)%tmpunique.length]].components
            var c = copiedvertices[tmpunique[(i + 2)%tmpunique.length]].components
            console.log("LOGGED POINT 1 POINT 2 POINT 3 AND LEFTOVER THING")
            if ((! arraysEqual(a, b)) & (! arraysEqual(a,c)) & (! arraysEqual(b,c))) {
                var constructed = this.constructtriangle(copiedvertices[tmpunique[i%tmpunique.length]], copiedvertices[tmpunique[(i + 1)%tmpunique.length]], copiedvertices[tmpunique[(i+2)%tmpunique.length]])
                firstvec = constructed[0]
                secondvec = constructed[1]
                thirdvec = constructed[2]
                midtocenter = constructed[3]
                center = constructed[4]
                if (this.checktriangle(firstvec, secondvec, midtocenter)) {
                    // this.garbage(center.components[0], center.components[1])
                    // this.garbage(constructed[5].components[0], constructed[5].components[1])
                    // this.garbage(constructed[7].components[0], constructed[7].components[1])
                    // this.garbage(constructed[6].components[0], constructed[6].components[1])
                    this.drawtriangle(constructed[5], constructed[6], constructed[7], center)
                    tmpunique.splice(((i+1)%tmpunique.length), 1) // remove the node hahahahaha let's go
                    console.log("SPLICING")
                } else {
                    i+=1
                }
            } else{
                i+=1
            }
        }
    }


    copy() { // this function should be used to get the screen components from the vertices
        var copied = this.vertices.map((value, index)=>{var tmp = new Vector(); tmp.components = JSON.parse(JSON.stringify(this.vertices[index].screencomponents)); return tmp;})
        console.log("COPYING ----------------------------------------------------------")
        console.log(copied)
        return copied
    } // need to have screen components, those depend on user moving, and need the vertices in proper order. This makes most sense.

    checktriangle(firstvec, secondvec, midtocenter) {
        console.log("CHECKING A TRIANGLE ___________________________________________")
        var a = Vector.getangle(firstvec, secondvec, this.orientation)
        console.log(a)
        var b = Vector.getangle(firstvec, midtocenter, this.orientation)
        console.log(b)
        console.log(firstvec)
        console.log(secondvec)
        console.log(midtocenter)
        // ADD CHECK TOO SEE IF THERE IS A POINT THAT IS CLOSE TO MID POINT THEN THE CENTER IN THAT REGION!!!!!
        if (b<=a) {
            console.log("YUPPER")
            return true
        } else if (b>a) {
            console.log("NOPER")
            return false
        }
    }

    checkorigins(base, height) {
        var basisvecs = [new Vector(1, 0), new Vector(0, 1)]
        var basedir = base.dot(basisvecs[0])
        var heightdir = height.dot(basisvecs[1])
        var change
        if ((basedir == 0) | (heightdir ==0)) {
            basedir = base.dot(basisvecs[1])
            heightdir = height.dot(basisvecs[0])
            base, height = height, base
            change = "Yes"
        } else {
            change = "No"
        }
        var origin
        var skeworigin
        var skewmult
        var bottommodpos
        var leftmodpos
        if (heightdir > 0) {
            origin = 'bottom'
            skeworigin = 'bottom'
            bottommodpos = 0
            skewmult = 1

        } else {
            origin = 'top'
            skeworigin = 'top'
            bottommodpos = 1
            skewmult = -1
        }
        if (basedir > 0) {
            origin = origin + ' left'
            skeworigin = skeworigin + ' right'
            skewmult = 1 * skewmult
            leftmodpos = 0
        } else {
            origin = origin + ' right'
            skeworigin = skeworigin + ' left'
            skewmult = -1 * skewmult
            leftmodpos = 1
        } // this should all work!

        // now calculate the angle, let's say the base matches with bassisvecs [0], then the angle is given by acos ( adj / hyp ) hyp = the base, adj is proj of base onto real base, which is given by 
        // proj(b', b) fake base ONTO real base = (b . b') x (b/|b|^2) = cos ß * |b| * |b'| * (b/|b|^2) = |proj| * b/|b| = proj
        var onto_base
        if (basedir>0) {
            onto_base = basisvecs[0]
        } else {
            onto_base = basisvecs[0].multiply(-1)
        }
        var adj =  (onto_base.multiply(1/Math.pow(onto_base.magnitude(),2))).multiply(base.dot(onto_base))
        var angle = Math.acos( adj.magnitude()/ base.magnitude() )
        angle = angle  * Math.sign(onto_base.components[0]) * (-Math.sign(base.dot(basisvecs[1])))
        return [origin, skeworigin, skewmult, bottommodpos, leftmodpos, angle, change]
    }


    drawrighttriangle(firstpoint, secondpoint, thirdpoint) { // third point is the middle point, need to get those coordinates to decide on the position
        var basevec = Vector.difference(thirdpoint, firstpoint)
        var heightvec = Vector.difference(thirdpoint, secondpoint)
        var bigdiv = document.createElement('div')
        var divorigin
        var skeworigin
        var skewmult
        var bottomomdpos 
        var leftmodpos
        
        var tmp = this.checkorigins(basevec, heightvec)
        divorigin = tmp[0]
        skeworigin = tmp[1]
        skewmult = tmp[2]
        bottomomdpos = tmp[3]
        leftmodpos = tmp[4]
        if (tmp[6] == "Yes") {
            basevec = Vector.difference(thirdpoint, secondpoint)
            heightvec = Vector.difference(thirdpoint, firstpoint)
        }
        var base = basevec.magnitude()
        var height = heightvec.magnitude()
        // console.log("BASE and height  TRIANGLE  ____________________________________________________________")
        // console.log(base)
        // console.log(height)
        // console.log("BASE  TRIANGLE  ____________________________________________________________")
        bigdiv.classList.add("projection")
        bigdiv.style.height= String(height) + 'px'
        bigdiv.style.width = String(base) + "px"
        bigdiv.style.overflow = 'hidden'
        var smalldiv = document.createElement('div')
        smalldiv.classList.add("projection")
        smalldiv.style.height= String(height) + 'px'
        smalldiv.style.width = String(base) + "px"
        bigdiv.style.transformOrigin = String(divorigin)
        smalldiv.style.transformOrigin = String(skeworigin)

        var angle = Math.atan(base/height)

        smalldiv.style.transform = "skewX(" + String(skewmult*angle) + "rad)" 

        smalldiv.style.backgroundColor = this.color
        var offset = tmp[5]
        bigdiv.style.rotate = String(offset) + "rad"
        bigdiv.style.position = 'absolute'
        bigdiv.appendChild(smalldiv)
        bigdiv.style.left = (thirdpoint.components[0] - base * leftmodpos) + 'px'
        bigdiv.style.bottom = (thirdpoint.components[1] - height * bottomomdpos) + 'px'
        mainframe.appendChild(bigdiv)
    }

    drawtriangle(firstpoint, secondpoint, thirdpoint, midpoint) {
        this.drawrighttriangle(firstpoint, secondpoint, midpoint)
        this.drawrighttriangle(thirdpoint, secondpoint, midpoint)
    }



    constructtriangle(firstpoint, secondpoint, thirdpoint){
        // if (Vector.difference(firstpoint, secondpoint).magnitude()<Vector.difference(firstpoint, thirdpoint).magnitude()) {
        //     secondpoint, thirdpoint = thirdpoint, secondpoint
        // }
        // var thirdvec = Vector.difference(thirdpoint, firstpoint)
        // var firstvec = Vector.difference(firstpoint, secondpoint)
        // var secondvec = Vector.difference(secondpoint, thirdpoint)
        // need to get center, so need to find base, could be any of the vecs, so need to get proj, check if is inside the triangle, if yes then valid base
        var points = [firstpoint, secondpoint, thirdpoint]
        for (let i=0; i<points.length; i++) {
            var vec = Vector.difference(points[(i % points.length)], points[((i+1) % points.length)])
            var onto = Vector.difference(points[(i % points.length)], points[((i+2) % points.length)])
            var projected = vec.proj(onto)
            if (((projected[1])<1) & (projected[1]>0)) {
                var firstvec = Vector.difference(points[i % points.length], points[(i+1) % points.length])
                var secondvec = Vector.difference(points[(i+1) % points.length], points[(i+2) % points.length])
                var thirdvec = Vector.difference(points[(i+2) % points.length], points[(i) % points.length])
                var firstpoint2 = points[i % points.length]
                var secondpoint2 = points[(i+1) % points.length]
                var thirdpoint2 = points[(i+2) % points.length]
                var tmp = firstvec.proj(Vector.difference(firstpoint2, thirdpoint2))

            }
            
        }
        // var tmp = firstvec.proj(Vector.difference(firstpoint, thirdpoint)) // should get middle of new triangle
        // var tmp = Vector.difference(firstpoint, secondpoint).multiply(0.5)
        var center =Vector.difference(firstpoint2, thirdpoint2).multiply(tmp[1])
        // var midtocenter = Vector.sum(midtocenter, secondpoint)
        var center = Vector.sum(center, firstpoint2)
        var midtocenter = Vector.difference(secondpoint2, center) // should get middle point to center vector


        return [firstvec, secondvec, thirdvec, midtocenter, center, firstpoint2, secondpoint2, thirdpoint2]
    } 

}

var world = new World(3)
class obsScreen {
    constructor(corners) {
        this.corners = corners // this will be relative to the observer so yeah
        this.basisx = null
        this.basisy = null
        this.basis = this.getbasis()
    }

    getbasis() {
        // the two basis vectors will be the bottom left [2] to bottom right [3 | -1] and bottom left [2] to top left [1]
        this.basisx = Vector.difference(this.corners[2], this.corners[3])
        this.basisy = Vector.difference(this.corners[2], this.corners[1])
        // this makes some sense
    }

    basischange(vector) { // takes in a point rel to screen, so vector from bottom left corner to point, and then converts it to basis
        var tmp = new Vector(vector.dot(this.basisx)/Math.pow(this.basisx.magnitude(),2), vector.dot(this.basisy)/Math.pow(this.basisy.magnitude(), 2))
        // console.log("CHECKING IF BASIS CHANGED TIMES BASIS IS EQUAL TO ORIGINAL VECTOR!")
        // console.log(vector)
        // console.log(tmp)
        // console.log(this.basisx.multiply(tmp.components[0]).add(this.basisy.multiply(tmp.components[1])))
        return tmp
        /* we do this as the dot product a.b = cos(ß) * ||a|| ||b||, so to get the ratio of this to the x basis, we have to divide by a^2 to get ||b||/||a|| * cos(ß)*/
    }
}
directiontracker = document.getElementById('direction')
class Observer {
    constructor (...coordinates) {
        this.position = new Vector(...coordinates)
        this.direction = new Vector(1,0,0)
        this.theta = 0
        this.phi = 0
        this.screendistance = 2
        this.thetarad = 30 // so 120 degress 
        this.phirad = 30 // so 100˚ of top down visual field isn't this value usually 60??? do they mean rad or total... anways
        this.corners = [] // this will store corners the the following order, top right, top left, bottom left, bottom right
        this.screen = this.constructscreen()
        console.log(this.screen.corners)
    }

    updatedirection(dphi=0, dtheta=0) {

        var newphi = (dphi/180 + this.phi) 
        var newtheta = (dtheta/180 + this.theta) 
        if ((dphi==0) & (dtheta==0)) {
            var direction = this.direction
        } else {
            var direction = new Vector()
            direction.components = JSON.parse(JSON.stringify(this.direction.components))// if I'm passing in params I want it to out a new vector so I can compute shit instead of updating the direction
            // direction.components = this.direction.components
        }
        direction.components[2] = Math.sin(newphi * Math.PI)
        direction.components[1] = Math.sin(newtheta * Math.PI) 
        direction.components[0] = Math.cos(newtheta * Math.PI) 
        var mag = direction.magnitude()
        direction.components.map((component, index) => {direction.components[index] = component})
        return direction
    }

    constructscreen() { //
        var m = Math.tan((this.thetarad/180) * Math.PI) * this.screendistance // perpendicular to x, y plane
        var n = Math.tan((this.phirad/180) * Math.PI) * this.screendistance
        var d = Math.sqrt((Math.pow(m,2) + Math.pow(n,2))) // get hyp
        var hyp = Math.sqrt(Math.pow(d,2) + Math.pow(this.screendistance, 2)) // this is the hyp magnitude to reach each corner and it is equal albeit in different directions
        var tmp = [[this.thetarad, this.phirad], [-this.thetarad, this.phirad], [-this.thetarad, -this.phirad], [this.thetarad, -this.phirad]]
        var points_rel_obs = tmp.map((corner, index) => {return this.updatedirection(corner[1], corner[0]).multiply(hyp)}) // returns coordinates of the screen corners relative to observer for ease, should perhaps add vertices, as how else will I do the dot product
        return new obsScreen(points_rel_obs) // this should build screen object with : corners, and basis vectors
        // var center = Vector.sum(this.position, this.direction.multiply(this.screendistance)) // Don't need center if I do I can always get it pretty easily from screen distance * direction
    }



    handleEvent(event) { // works wow easy ok not so easy lol
        var previousprojections = document.getElementsByClassName('projection')
        for (let i=0; i<previousprojections.length; i++) {
            previousprojections[0].remove()
        }
        if (event.type == 'mousemove') {
            this.mousemovehandling(event)
        } else if (event.type == 'keydown') {
            this.keyhandling(event)
        }
    }

    mousemovehandling (event) {
        this.theta = (event.movementX/screen.availWidth + this.theta) % 2
        if (((-event.movementY/screen.availWidth + this.phi) % 2)>1/2) {
            this.phi = 1/2
        } else if (((-event.movementY/screen.availWidth + this.phi) % 2)<(-1/2)) {
            this.phi = -1/2
        } else {
            this.phi = (-event.movementY/screen.availWidth + this.phi) % 2
        }
        this.updatedirection()
        directiontracker.innerHTML = this.direction.components
        this.screen = this.constructscreen()
        world.project(me)
    }

    async keyhandling(event) {
        if (event.key == 'w') {
            for (let i=0; i<20; i++) {
                me.position.add(this.direction.multiply(0.1))
                me.position.components[2] = 0
                for (let j=0; j<me.position.components.length; j++) {
                    positiontracker.childNodes[(2*j + 1)].innerHTML = positiontracker.childNodes[(2*j + 1)].innerHTML.slice(0, 14) + String(Number.parseFloat(me.position.components[j]).toFixed(1))
                }
                this.updatedirection()
                directiontracker.innerHTML = this.direction.components
                this.screen = this.constructscreen()
                world.project(me)
                await new Promise(r => setTimeout(r, 10));
            }
        }
        if (event.key == 's') {
            for (let i=0; i<20; i++) {
                me.position.add(this.direction.multiply(-0.1))
                me.position.components[2] = 0
                for (let j=0; j<me.position.components.length; j++) {
                    positiontracker.childNodes[(2*j + 1)].innerHTML = positiontracker.childNodes[(2*j + 1)].innerHTML.slice(0, 14) + String(Number.parseFloat(me.position.components[j]).toFixed(1))
                }
                this.updatedirection()
                this.screen = this.constructscreen()
                world.project(me)
                await new Promise(r => setTimeout(r, 10));
            }
        }
        if (event.key == 'a') {
             var tmp = new Vector()
            tmp.components = JSON.parse(JSON.stringify(this.direction.components))
            tmp.components[2] = 0
            var temporary =  tmp.components[1]
            var temporary2 =  -1* tmp.components[0]
            tmp.components[1] = temporary2
            tmp.components[0] = temporary
            for (let i=0; i<20; i++) {
                me.position.add(tmp.multiply(0.1))
                for (let j=0; j<me.position.components.length; j++) {
                    positiontracker.childNodes[(2*j + 1)].innerHTML = positiontracker.childNodes[(2*j + 1)].innerHTML.slice(0, 14) + String(Number.parseFloat(me.position.components[j]).toFixed(1))
                }
                this.updatedirection()
                this.screen = this.constructscreen()
                world.project(me)
                await new Promise(r => setTimeout(r, 10));
            }
        }
        if (event.key == 'd') {
            var tmp = new Vector()
            tmp.components = JSON.parse(JSON.stringify(this.direction.components))
            tmp.components[2] = 0
            var temporary =   -1*tmp.components[1]
            var temporary2 =   tmp.components[0]
            tmp.components[1] = temporary2
            tmp.components[0] = temporary
            for (let i=0; i<20; i++) {
                me.position.add(tmp.multiply(0.1))
                for (let j=0; j<me.position.components.length; j++) {
                    positiontracker.childNodes[(2*j + 1)].innerHTML = positiontracker.childNodes[(2*j + 1)].innerHTML.slice(0, 14) + String(Number.parseFloat(me.position.components[j]).toFixed(1))
                }
                this.updatedirection()
                this.screen = this.constructscreen()
                world.project(me)
                await new Promise(r => setTimeout(r, 10));
            }
        }
    }

}

function cloneVertex(vertex) {
    // return new Vertex(JSON.parse(JSON.stringify(vertex.components)))    
    return vertex
}



const a = new Vector(1,0)
friction_const = 0.1 
friction_coef = 0.7
epsilon = 1e-1



var vertex1 = new Vertex(0, 0, 0)
var vertex2 = new Vertex(0,0,1)
var vertex3 = new Vertex(0,1,0)
var vertex4 = new Vertex(1,0,0)
var vertex5 = new Vertex(0,1,1)
var vertex6 = new Vertex(1,0,1)
var vertex7 = new Vertex(1,1,0)
var vertex8 = new Vertex(1,1,1)

vertex3.addneighbour(vertex1, vertex7, vertex5)
vertex4.addneighbour(vertex1, vertex7, vertex6)
vertex5.addneighbour(vertex8, vertex3, vertex2)
vertex6.addneighbour(vertex4, vertex2, vertex8)
vertex7.addneighbour(vertex3, vertex4, vertex8)
vertex8.addneighbour(vertex7, vertex5, vertex6)
vertex2.addneighbour(vertex1, vertex5, vertex6)
vertex1.addneighbour(vertex2, vertex3, vertex4)


// var vertex1 = new Vertex(0, 0, 0)
// var vertex2 = new Vertex(0,0,1)
// var vertex3 = new Vertex(0,1,0)
// var vertex4 = new Vertex(0,1,1)

// vertex3.addneighbour(vertex1, vertex4)
// vertex2.addneighbour(vertex4, vertex1)
// vertex1.addneighbour(vertex2, vertex3)
// vertex4.addneighbour(vertex3, vertex2)
// // var face1 = new face('blue', new Vertex(0,0), new Vertex(10, 0), new Vertex (10,10 ),new Vertex(0,10))
// // var face2 = new face('blue', new Vertex(0,0), new Vertex(0,10), new Vertex (10,10 ), new Vertex(10, 0))

// var face1 = new face('red', cloneVertex(vertex1), cloneVertex(vertex2), cloneVertex(vertex4), cloneVertex(vertex3))

// var vertex1 = new Vertex(0, 0, 0)
// var vertex2 = new Vertex(0,0,1)
// var vertex3 = new Vertex(0,1,0)
// var vertex4 = new Vertex(0,1,1)

// vertex3.addneighbour(vertex1, vertex4)
// vertex4.addneighbour(vertex2, vertex3)
// vertex2.addneighbour(vertex1, vertex2)
// vertex1.addneighbour(vertex2, vertex3)
// This should form a cube

// var cube = new object(vertex1, vertex3, vertex2, vertex4) // yay cube, maybe should add checkers but lowkey need results to keep working so need to seeee it weork. this is cool
var cube = new object(vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8)

var face2 = new face('red', cloneVertex(vertex1), cloneVertex(vertex2), cloneVertex(vertex6), cloneVertex(vertex4))
var face3 = new face('yellow', cloneVertex(vertex1), cloneVertex(vertex3), cloneVertex(vertex7), cloneVertex(vertex4))
var face4 = new face('orange', cloneVertex(vertex8), cloneVertex(vertex6), cloneVertex(vertex2), cloneVertex(vertex5))
var face5 = new face('purple', cloneVertex(vertex7), cloneVertex(vertex8), cloneVertex(vertex5), cloneVertex(vertex3))
var face6 = new face('green', cloneVertex(vertex7), cloneVertex(vertex8), cloneVertex(vertex6), cloneVertex(vertex4))
cube.faces = [face2, face3, face4, face5, face6]

// object.faces = [face1, face2, face3, face4, face5, face6]
// cube.faces = [face1]
// var cube = new object(vertex1, vertex2, vertex3, vertex4) // yay cube, maybe should add checkers but lowkey need results to keep working so need to seeee it weork. this is cool
world.addobject(cube)


function friction(x) {
    if (x>friction_const) {
        return ((x - friction_const) * friction_coef)
    } else if (x>epsilon) {
        x = 0
    }
}

mainframe.addEventListener('click', async (event) => {promise = (mainframe.requestPointerLock()).then(addEventListener('mousemove', me))}) // works surprisingly
var me = new Observer(5,0,0)

var positiontracker = document.getElementById('position')

addEventListener('keydown', me)

// addEventListener('mousemove', me) // works
// ArrowLeft ArrowUp ArrowDown ArrowRight
// w a s d
// box.style.left = String(event.pageX -parseInt( box.style.width)/2)+ 'px'; box.style.top= String(event.pageY - parseInt(box.style.height)/2 )+ 'px';


// Shape : closed shape has vertices, that link together, don't overlap and whatnot
// square : extends shape, must have 4 vertices, any two adjancent vectors must be perpendicular
// observer : rectangle, cone, position, direction, theta = angle in x,y phi = angle in x, z plane
// volume : a bunch of shapes binding together